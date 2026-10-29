"""A SeedLink 3 server for the tests, serving the records of miniSEED files.

Serves the 512-byte records of the given files the way a datalogger releases
them, in the order of the time of their last sample, at a chosen speed (the
first END starts the clock), to one client after another. It speaks the subset
of SeedLink 3 that ayzek's client uses, in multi-station mode: HELLO, STATION
sta [net], SELECT pattern, DATA [seq], END, INFO and BYE. Sequence numbers are
per station, in delivery order, 24 bit and wrapping; DATA with a number
resumes from that packet.

Faults, for tests/seedlink_check.sh:
  --drop-after N       drop the connection after N packets, halfway through
                       the next one, with a reset (repeatable: the Nth value
                       applies to the Nth connection)
  --stall-after N:S[:mute]
                       send nothing for S seconds after N packets in all, the
                       clock running on, then catch up; with :mute, ignore
                       INFO meanwhile (repeatable)
  --mute-info          ignore INFO (a client's keepalive gets no reply)
  --refuse-info        answer INFO with ERROR
  --info-packets N     split each INFO reply over N packets ("SLINFO *")
  --resume-overlap N   on DATA with a number, send the N packets before it
                       again, as a server that resumes inclusively would
  --seq-start HEX      every station's first sequence number (FFFFF0 tests
                       the wrap)
  --delay STA=S        station STA's records become available S seconds late
  --reorder N          every Nth record of each channel arrives after the
                       next one of its channel
  --exit-when-done     close and exit once a client has had every record

Standard library only.

    python3 tools/seedlink_mock.py --port 0 --port-file PORT --speed 20 \\
        --exit-when-done data/demo_ko/*.mseed
"""
import argparse
import datetime as dt
import os
import select
import socket
import struct
import sys
import threading
import time

REC = 512
MASK = 0xFFFFFF


def log(msg):
    print(f"mock {time.strftime('%H:%M:%S')} {msg}", file=sys.stderr, flush=True)


def sample_rate(factor, mult):
    # As src/mseed.hpp, so that end times order as ayzek's replay orders them.
    if factor == 0 or mult == 0:
        return 0.0
    if factor > 0 and mult > 0:
        return float(factor * mult)
    if factor > 0 and mult < 0:
        return -factor / mult
    if factor < 0 and mult > 0:
        return -mult / factor
    return 1.0 / (factor * mult)


class Record:
    __slots__ = ("raw", "net", "sta", "loc", "cha", "end", "avail", "seq")

    def __init__(self, raw):
        self.raw = raw
        self.sta = raw[8:13].decode("ascii", "replace").strip()
        self.loc = raw[13:15].decode("ascii", "replace").strip()
        self.cha = raw[15:18].decode("ascii", "replace").strip()
        self.net = raw[18:20].decode("ascii", "replace").strip()
        year, doy, hh, mm, ss, _, fract, n, factor, mult = struct.unpack(
            ">HHBBBBHHhh", raw[20:36])
        activity, = struct.unpack(">B", raw[36:37])
        correction, = struct.unpack(">i", raw[40:44])
        micro = 0
        off, = struct.unpack(">H", raw[46:48])
        for _ in range(16):
            if off == 0 or off + 8 > REC:
                break
            typ, nxt = struct.unpack(">HH", raw[off:off + 4])
            if typ == 1001:
                micro = struct.unpack(">b", raw[off + 5:off + 6])[0]
            off = nxt
        day = dt.date(year, 1, 1).toordinal() + doy - 1 - dt.date(1970, 1, 1).toordinal()
        ns = (((day * 24 + hh) * 60 + mm) * 60 + ss) * 10**9 + fract * 100_000 + micro * 1000
        if not activity & 0x02:
            ns += correction * 100_000
        rate = sample_rate(factor, mult)
        self.end = ns / 1e9 + (n / rate if rate else 0.0)
        self.avail = self.end
        self.seq = 0


def load(paths, delays, reorder, seq_start):
    records = []
    for p in paths:
        data = open(p, "rb").read()
        for i in range(0, len(data) - REC + 1, REC):
            records.append(Record(data[i:i + REC]))
    if not records:
        sys.exit("seedlink_mock: no records")
    for r in records:
        r.avail = r.end + delays.get(r.sta, 0.0)
    records.sort(key=lambda r: r.avail)  # stable: file order breaks ties
    if reorder:
        count, i = {}, 0
        while i < len(records) - 1:
            r = records[i]
            key = (r.net, r.sta, r.loc, r.cha)
            count[key] = count.get(key, 0) + 1
            if count[key] % reorder == 0:
                j = next((j for j in range(i + 1, len(records))
                          if (records[j].net, records[j].sta, records[j].loc,
                              records[j].cha) == key), None)
                if j is not None:  # r goes just after the next of its channel
                    records.insert(j, records.pop(i))
                    continue
            i += 1
    by_station = {}
    for r in records:
        lst = by_station.setdefault((r.net, r.sta), [])
        r.seq = (seq_start + len(lst)) & MASK
        lst.append(r)
    return records, by_station


def selector_matches(pattern, rec):
    """True/False if the pattern selects or excludes the record, None if it
    does not concern it (a negative pattern that does not match)."""
    neg = pattern.startswith("!")
    p = pattern[1:] if neg else pattern
    p, _, kind = p.partition(".")
    if kind and kind.upper() != "D":
        return None if neg else False
    loc, cha = (None, p) if len(p) <= 3 else (p[:-3], p[-3:])

    def like(pat, s):
        return len(pat) == len(s) and all(a == "?" or a == b for a, b in zip(pat, s))
    hit = like(cha, rec.cha) and (loc is None or like(loc.replace("-", " "), rec.loc.ljust(2)))
    if neg:
        return False if hit else None
    return hit


def selected(selectors, rec):
    if not selectors:
        return True
    pos = [s for s in selectors if not s.startswith("!")]
    if any(selector_matches(s, rec) is False for s in selectors if s.startswith("!")):
        return False
    return not pos or any(selector_matches(s, rec) for s in pos)


def info_packets(xml, n):
    """The INFO reply as `n` SLINFO packets, each a miniSEED log record."""
    chunk = max(1, -(-len(xml) // n))
    parts = [xml[i:i + chunk] for i in range(0, len(xml), chunk)][:n]
    now = dt.datetime.now(dt.UTC)
    out = []
    for i, part in enumerate(parts):
        hdr = bytearray(REC)
        hdr[0:8] = b"000000D "
        hdr[8:20] = b"INFO " + b"  " + b"LOG" + b"SL"  # station, location, channel, network
        struct.pack_into(">HHBBBBH", hdr, 20, now.year, now.timetuple().tm_yday,
                         now.hour, now.minute, now.second, 0, now.microsecond // 100)
        struct.pack_into(">HhhBBBBihH", hdr, 30, len(part), 0, 0, 0, 0, 0, 1, 0, 64, 48)
        struct.pack_into(">HHBBBB", hdr, 48, 1000, 0, 0, 1, 9, 0)
        hdr[64:64 + len(part)] = part
        out.append((b"SLINFO *" if i < len(parts) - 1 else b"SLINFO  ") + bytes(hdr))
    return out


class Server:
    def __init__(self, a, records, by_station):
        self.a = a
        self.records = records
        self.by_station = by_station
        self.lock = threading.Lock()
        self.wall0 = None
        self.t0 = records[0].avail
        self.sent_total = 0
        self.stalls = sorted(a.stall)  # (packets, seconds, mute), each once
        self.connections = 0

    def now(self):
        """Stream time: records with avail <= now() may be sent."""
        if self.a.speed <= 0:
            return float("inf")
        with self.lock:
            if self.wall0 is None:
                return self.t0
            return self.t0 + (time.monotonic() - self.wall0) * self.a.speed

    def start_clock(self):
        with self.lock:
            if self.wall0 is None:
                self.wall0 = time.monotonic()
                log("clock started")

    def serve(self, conn, number):
        try:
            Session(self, conn, number).run()
        except (ConnectionError, OSError) as e:
            log(f"connection {number}: {e}")
        finally:
            try:
                conn.close()
            except OSError:
                pass


class Session:
    def __init__(self, server, conn, number):
        self.s = server
        self.a = server.a
        self.conn = conn
        self.number = number
        self.buf = b""
        drops = self.a.drop_after
        self.drop_after = drops[number - 1] if number <= len(drops) else None
        self.muted = False

    def send(self, data):
        self.conn.sendall(data)

    def read_line(self, timeout=None):
        while True:
            for sep in (b"\r\n", b"\n", b"\r"):
                i = self.buf.find(sep)
                if i >= 0:
                    line, self.buf = self.buf[:i], self.buf[i + len(sep):]
                    return line.decode("ascii", "replace").strip()
            if timeout is not None:
                r, _, _ = select.select([self.conn], [], [], timeout)
                if not r:
                    return None
            data = self.conn.recv(4096)
            if not data:
                raise ConnectionError("closed by the client")
            self.buf += data

    def info(self, args):
        if self.a.mute_info or self.muted:
            return
        if self.a.refuse_info:
            self.send(b"ERROR\r\n")
            return
        level = args[0].upper() if args else "ID"
        xml = (f'<?xml version="1.0"?>\n<seedlink software="SeedLink v3.1 (ayzek mock)" '
               f'organization="ayzek test server" level="{level}"/>\n').encode()
        for p in info_packets(xml, self.a.info_packets):
            self.send(p)

    def run(self):
        stations, current = {}, None
        while True:  # handshake
            line = self.read_line()
            if not line:
                continue
            cmd, *args = line.split()
            cmd = cmd.upper()
            if cmd == "HELLO":
                self.send(b"SeedLink v3.1 (ayzek mock) :: SLPROTO:3.1\r\nayzek test server\r\n")
            elif cmd == "STATION" and args:
                key = next((k for k in self.s.by_station
                            if k[1] == args[0].upper() and (len(args) < 2 or k[0] == args[1].upper())),
                           None)
                current = key
                if key is None:
                    self.send(b"ERROR\r\n")
                else:
                    stations[key] = {"selectors": [], "start": None}
                    self.send(b"OK\r\n")
            elif cmd == "SELECT":
                if current is None or len(args) > 1 or (args and not all(
                        c.isalnum() or c in "?!.-" for c in args[0])):
                    self.send(b"ERROR\r\n")
                else:
                    if args:
                        stations[current]["selectors"].append(args[0].upper())
                    else:
                        stations[current]["selectors"].clear()
                    self.send(b"OK\r\n")
            elif cmd == "DATA":
                if current is None:
                    self.send(b"ERROR\r\n")
                else:
                    stations[current]["start"] = int(args[0], 16) if args else None
                    self.send(b"OK\r\n")
            elif cmd == "END":
                break
            elif cmd == "INFO":
                self.info(args)
            elif cmd == "BYE":
                return
            else:
                self.send(b"ERROR\r\n")
        if not stations:
            return
        log(f"connection {self.number}: " + ", ".join(
            f"{k[0]}.{k[1]} from {'%06X' % v['start'] if v['start'] is not None else 'now'}"
            for k, v in stations.items()))
        now = self.s.now()  # before the first END starts the clock: its start
        self.s.start_clock()
        self.stream(stations, now)

    def stream(self, stations, now):
        # Index within each station's records of the first one to send: the
        # one numbered as asked, or the first not yet available.
        for key, st in stations.items():
            recs = self.s.by_station[key]
            if st["start"] is not None:
                k = (st["start"] - recs[0].seq) & MASK
                k = min(k, len(recs))
                if self.a.resume_overlap:
                    k = max(0, k - self.a.resume_overlap)
            else:
                k = next((i for i, r in enumerate(recs) if r.avail >= now), len(recs))
            st["k"] = k
        index = {key: 0 for key in stations}
        sent = 0
        for rec in self.s.records:
            key = (rec.net, rec.sta)
            if key not in stations:
                continue
            k = index[key]
            index[key] += 1
            st = stations[key]
            if k < st["k"] or not selected(st["selectors"], rec):
                continue
            self.wait_for(rec)
            packet = b"SL%06X" % rec.seq + rec.raw
            if self.drop_after is not None and sent >= self.drop_after:
                self.send(packet[:260])
                log(f"connection {self.number}: dropped after {sent} packets")
                self.conn.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER,
                                     struct.pack("ii", 1, 0))
                return
            self.send(packet)
            sent += 1
            self.after_packet()
        log(f"connection {self.number}: every record sent ({sent} packets)")
        if self.a.exit_when_done:
            self.conn.shutdown(socket.SHUT_WR)
            time.sleep(0.5)
            self.conn.close()
            log("exiting")
            os._exit(0)
        while True:  # idle: answer INFO until the client goes
            self.poll_commands(1.0)

    def after_packet(self):
        with self.s.lock:
            self.s.sent_total += 1
            stall = None
            if self.s.stalls and self.s.sent_total >= self.s.stalls[0][0]:
                stall = self.s.stalls.pop(0)
        if stall:
            _, seconds, mute = stall
            log(f"stalling {seconds:g} s after {self.s.sent_total} packets"
                f"{', ignoring INFO' if mute else ''}")
            self.muted = mute
            end = time.monotonic() + seconds
            while time.monotonic() < end:
                self.poll_commands(min(0.1, end - time.monotonic()))
            self.muted = False
            log("stall over")

    def wait_for(self, rec):
        while True:
            now = self.s.now()
            if rec.avail <= now:
                return
            self.poll_commands(min(0.05, (rec.avail - now) / self.a.speed))

    def poll_commands(self, timeout):
        line = self.read_line(max(timeout, 0.0))
        if line is None:
            return
        cmd, *args = line.split() or [""]
        if cmd.upper() == "INFO":
            self.info(args)
        elif cmd.upper() == "BYE":
            raise ConnectionError("BYE")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=0, help="0: any free port")
    ap.add_argument("--port-file", help="write the port here once listening")
    ap.add_argument("--speed", type=float, default=1.0, help="0: as fast as possible")
    ap.add_argument("--drop-after", type=int, action="append", default=[])
    ap.add_argument("--stall-after", dest="stall", action="append", default=[])
    ap.add_argument("--mute-info", action="store_true")
    ap.add_argument("--refuse-info", action="store_true")
    ap.add_argument("--info-packets", type=int, default=1)
    ap.add_argument("--resume-overlap", type=int, default=0)
    ap.add_argument("--seq-start", default="0")
    ap.add_argument("--delay", action="append", default=[], help="STA=SECONDS")
    ap.add_argument("--reorder", type=int, default=0)
    ap.add_argument("--exit-when-done", action="store_true")
    a = ap.parse_args()
    stalls = []
    for spec in a.stall:
        n, secs, *flag = spec.split(":")
        stalls.append((int(n), float(secs), flag == ["mute"]))
    a.stall = stalls
    delays = {d.split("=")[0].upper(): float(d.split("=")[1]) for d in a.delay}
    records, by_station = load(a.files, delays, a.reorder, int(a.seq_start, 16))
    log(f"{len(records)} records of {len(by_station)} stations, "
        f"{records[0].avail:.2f} to {records[-1].avail:.2f}, speed {a.speed:g}")

    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    srv.bind((a.host, a.port))
    srv.listen(4)
    port = srv.getsockname()[1]
    if a.port_file:
        tmp = a.port_file + ".tmp"
        with open(tmp, "w") as f:
            f.write(f"{port}\n")
        os.replace(tmp, a.port_file)
    log(f"listening on {a.host}:{port}")
    server = Server(a, records, by_station)
    while True:
        conn, peer = srv.accept()
        server.connections += 1
        log(f"connection {server.connections} from {peer[0]}:{peer[1]}")
        threading.Thread(target=server.serve, args=(conn, server.connections),
                         daemon=True).start()


if __name__ == "__main__":
    main()

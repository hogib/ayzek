# 19 · Live input: SeedLink

`--seedlink HOST:PORT --streams NET.STA[,...]` takes the place of the
miniSEED files: one thread receives every station's records from a SeedLink
server and releases each into its station's rings as it arrives. Everything
after the rings — processors, network stage, report — is the code a replay
runs, and so are decoding, reordering and the ring writes, which both paths
reach through `StationFeed` (`src/pipeline/ingest.hpp`).

```
SeedLink server ── TCP ──> Client (src/seedlink.cpp)            one thread
                             HELLO, STATION/SELECT/DATA per station, END
                             packets -> duplicate check -> per station
                           LiveInput (src/pipeline/live.cpp)
                             archive · checks · drop if a ring behind
                             StationFeed: decode -> Reorderer(0) -> rings ──> processors (as replay)
```

| file | |
|---|---|
| `src/seedlink.hpp`, `src/seedlink.cpp` | the protocol: packet framing, handshake, sequence numbers, the client with its timeouts and reconnection |
| `src/pipeline/live.hpp`, `src/pipeline/live.cpp` | records to stations: what is accepted, what is rejected and why, the archive, drops |
| `src/pipeline/ingest.hpp` | `StationFeed`, shared with the replay |
| `app/ayzek.cpp` | options, the network stage's ordering for live input, the report |
| `tools/seedlink_mock.py` | a SeedLink server for the tests, standard library only |

## Protocol

SeedLink 3, written from the SeisComP description of the server
(seiscomp.de, `apps/seedlink`) and, for what version 4 keeps of version 3's
semantics, the FDSN SeedLink 4 specification. No SeedLink library is used.
The client uses multi-station mode only:

```
HELLO                  two lines: server version, organisation
STATION sta net        OK or ERROR     ┐
SELECT HH?.D           OK or ERROR     ├ per station
DATA [seq]             OK or ERROR     ┘
END                    no reply; packets follow
```

- **Packets.** "SL", six hexadecimal digits of sequence number, one 512-byte
  miniSEED record: 520 bytes. `FrameReader` splits the byte stream into
  packets however TCP cuts it, and anything that is not a packet ends the
  connection, since a stream cannot be resynchronised after garbage.
- **INFO.** The only command allowed once data flows. Replies are packets
  headed `SLINFO` whose 8th byte is `*` while more of the reply follows and a
  space on the last; a refused INFO is answered `ERROR\r\n` in the packet
  stream. The client sends `INFO ID` as its keepalive and discards the reply.
  `END` (the end of a dial-up transfer, which DATA never asks for) is read as
  the server ending the session.
- **Selectors.** `HH?.D`: the HH channels' data records, any location, or
  `LLHH?.D` for `NET.STA.LL`. The empty location (`NET.STA.--`) has no
  selector spelling every server reads alike, so it is selected client-side.
  Other channels and rates are filtered by the client after that (below).
- **Sequence numbers** are per station, 24 bit, wrapping after FFFFFF.
  `DATA n` asks for the packets from n on, so a station resumes after its last
  packet s with `DATA s+1`; `DATA` alone starts with the next packet, real
  time, which is what the first connection asks for. With selectors the
  numbers have gaps, so a gap in them says nothing about loss.
- **Checked against a real server** (GEOFON, `geofon.gfz.de:18000`, station
  GE.WLF, 2026-10-29): it is a SeedLink 4 server serving v3 clients (HELLO
  answers `SeedLink v4.0 [HMB SeedLink v0.2 (2026.110)] :: SLPROTO:4.0`);
  the INFO ID reply is one packet headed `SLINFO  `; and `DATA n` sends
  packet n again, so `last+1` is the right resume. A 90 s run through a local
  proxy that reset the connection after 30 s reconnected after 1 s, and the
  archive's 21 records per channel have no gap or overlap and no duplicate
  across the reconnection; latency (arrival minus the record's end) had a
  median of 0.9 s.

## Failure handling

- **Timeouts.** Connecting: 10 s per address. Every handshake reply:
  `--seedlink-timeout` (120 s). Data: when nothing at all has been received
  for `--seedlink-keepalive` (60 s), `INFO ID` is sent; its reply keeps a
  quiet but healthy connection alive, and when nothing has been received for
  `--seedlink-timeout` the connection is dropped. With continuous 100 Hz data
  no keepalive is ever sent. IRIS asks for no more than one keepalive every 4
  minutes on an idle connection (`--seedlink-keepalive 240
  --seedlink-timeout 300` there).
- **Reconnection.** After any failure — connection refused or reset, a
  timeout, a protocol error, the server closing — the client reconnects after
  1, 2, 4, ... up to 30 s, back to 1 s after a session that delivered data,
  and every station resumes after its last packet, so nothing the server
  still holds is lost. A station that has not delivered a packet yet has no
  number to resume from and asks again for real-time data, so what the
  server buffered for it during the outage is not received.
- **Duplicates.** A packet whose number is not after the station's last is
  discarded as already delivered, unless its data is newer than anything of
  its channel delivered: then the server's numbering has started again (a
  restarted server) and the client follows the new numbering, with a warning.
  So a server that resumes inclusively, or sends a packet twice, does not
  duplicate data, and one that restarts its numbering does not stop the
  station.
- **Refusals.** ERROR to a station's STATION, SELECT or DATA in the first
  session ends the run with that message: the stream list is wrong. In a
  later session the station is left out until the next connection, with a
  warning.
- **Records.** The pipeline takes HHZ, HHN and HHE at exactly 100 Hz. Every
  other record is counted by reason and reported, the first of each reason as
  it happens: another channel; HH1/HH2, which would need a rotation live
  input does not do; another rate (there is no resampling); a second location
  code of an any-location station (that of the first HHZ, HHN or HHE record
  at 100 Hz is kept, `NET.STA.LOC` chooses); an encoding the decoder does not support. A record
  that does not decode is counted too (`bad data`).
- **A full ring.** The reader never waits for a processor. While a
  component holds more than a ring's worth of samples back
  (`kMaxQueuedSamples`), its records are dropped, counted and reported. A
  dropped record becomes a gap when the next one arrives (the reorderer
  declares it), so no later sample moves in time; a replay instead holds the
  source back, which a live server cannot be. The limit is per component, not
  per station: during an outage of one component the processor cannot
  advance, the other two fill their rings, and when the missing one returns
  it holds nothing back, so its records are taken and its outage becomes a
  gap. Counted per station, its records would be dropped too and the station
  never recover.
- **Stopping.** The first SIGINT or SIGTERM, or the end of `--duration`,
  stops the input; what has been received is processed and the normal report
  printed, followed by the input's section. A second signal stops at once.
  An error after the input has started (an output that cannot be written,
  say) stops the input and the processors and ends the run with that error.

## The network stage, live

In a replay the network stage takes messages in stream-time order, and only up
to the oldest promise of all stations (`05-pipeline.md`, *Network*), which
makes the result independent of thread timing. Live input keeps that rule,
which costs the alarm the lag of the slowest station: a message waits until
every station has processed its data past the message's time. That is a
record's duration (about 2 to 7 s for 100 Hz Steim2 at the KO stations) plus
telemetry delay. A station that stops sending would hold every alarm back
forever, so a message waits about `--live-wait` seconds (10); the stations
that held it back that long are then not waited for until they have caught
up with what has been processed, each time with a warning. What such a
station sends for earlier times is processed out of order, and the report
counts these messages. With none out of order the network stage saw exactly
the replay's order. The wait is timed from when the oldest held message
arrived, so a message that arrives later for an earlier time starts it again:
a message can wait longer than `--live-wait`, but not for ever, and none is
released out of order without being counted.

`ISSUED` lines give each alarm's wall-clock time and how long after the data
it came, which is what a live run measures and a replay cannot.

## Archive

`--archive DIR` appends every record the client delivers, as received and
before any check, to `DIR/STATION.mseed`, and flushes it after each record.
`ayzek --speed 0 DIR/*.mseed` replays a session. The replay gives the live
run's alarms when the live run dropped no record, processed nothing out of
order, and received no record late: the replay orders records by end time, as
a datalogger sends them, while live input takes them as they come, and a
record that arrives after a later one of its channel is discarded as stale
there (`05-pipeline.md`, *Reordering*). Rejected records are in the archive
too; a replay skips them for the same reasons, except HH1/HH2, which the
replay reads as N and E.

## Decoding

Live data comes from other networks' dataloggers, so the decoder gained
Steim1 (encoding 10) beside Steim2 and the 16- and 32-bit integers it had.
`tests/test_mseed.cpp` decodes one small file of each encoding, and of each
integer byte order, written by ObsPy (`tools/make_mseed_fixtures.py`, its
provenance in `tests/reference/mseed/fixtures.json`), and requires ObsPy's
own reading sample for sample: Steim1 with 8-, 16- and 32-bit differences,
Steim2 with 4- to 30-bit ones, integers within ±2^30. The little-endian
integer files have a big-endian header, as CATL's records do; records whose
header is little-endian too are not supported, and are counted as such.

## Tests

| test | what it checks |
|---|---|
| `test_seedlink` | packet framing at every split, INFO and ERROR frames, garbage; replies and commands; stream lists and addresses; sequence numbers across the wrap; the duplicate and restart rule; every rejection reason; the drop of records a ring behind and the gap they leave; the archive |
| `test_steim` | Steim1's three codes, its integrity check, the encodings and byte orders the header parser accepts |
| `test_mseed` | every encoding against ObsPy (above) |
| `seedlink_live` | ayzek against the mock serving `data/demo_ko` at 20x real time: the report and `--record` file equal the file replay's, and so do those of a replay of the archive; nothing dropped or out of order |
| `test_live_outage` | an outage of one component (HHE) longer than the ring: when it returns, its records are taken and the station recovers |
| `live_error_exit` | an error after the threads have started (an unwritable `--record`, against a closed local port) ends the live run with the error |
| `seedlink_reconnect` | the same with three dropped connections (halfway through a packet, with a reset), a stall answered by two-packet INFO replies to keepalives, a stall with INFO ignored past the read timeout, packets sent again after every resume and numbers across the wrap: the archive holds every record of the files exactly once, and the report and `--record` file equal the file replay's |

The live tests need `data/demo_ko` and python3 and skip without them; they run
natively only. `AYZEK_LIVE_SPEED` changes the speed.

**Blinding.** The held-out tests use live and recent data of KO, HL, HT, HA,
HP, HC and NZ, so no test or smoke test connects to KOERI, NOA or GeoNet, or
asks any server for those networks. Everything is tested against the mock on
development files.

## Not supported

- Rotation of HH1/HH2 to N/E, other sample rates, other channels.
- SeedLink 4 (`SLPROTO 4.0`, which a v4 server needs only for its new
  features), BATCH mode, uni-station mode, FETCH (dial-up) and TIME (time
  windows), INFO levels other than ID.
- Starting from the server's buffer: the first connection takes data from the
  next packet on, and sequence numbers are not saved across runs.
- `--assess`, which reads the replay files after the run (replay an archive
  instead), and `--speed`.
- Windows: the client needs POSIX sockets; `--seedlink` says so there.
- A catalogue (`--catalog`) is read, but only for testing: the report keeps
  the events within the span received, as a replay would.

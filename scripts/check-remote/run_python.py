import json
import os
import struct
import sys
from pathlib import Path
from urllib.request import urlopen

MASS = 90.05550
FROM_RT = 0.0
TO_RT = 5.0
EIC_PPM = 20.0
EIC_MZ = 0.005


def same_bits(left, right) -> bool:
    return len(left) == len(right) and all(
        struct.pack(">d", a) == struct.pack(">d", b) for a, b in zip(left, right)
    )


def stats(url: str, path: str) -> dict:
    origin = "/".join(url.split("/")[:3])
    with urlopen(origin + path) as response:
        return json.loads(response.read())


def main() -> int:
    url = os.environ.get("QUANTION_REMOTE_URL", "")
    if not url:
        print("QUANTION_REMOTE_URL is not set", file=sys.stderr)
        return 1

    root = Path(__file__).resolve().parents[2]
    fixture = root / "core/tests/fixtures/api/api.ion"
    sys.path.insert(0, str(root / "wrappers/python"))
    import quantion

    stats(url, "/reset")
    remote = quantion.parse_ion_remote(url)
    opening = stats(url, "/stats")

    stats(url, "/reset")
    eic_remote = quantion.calculate_eic(remote, MASS, FROM_RT, TO_RT, EIC_PPM, EIC_MZ)
    one_query = stats(url, "/stats")

    local = quantion.parse_ion_path(str(fixture))
    eic_local = quantion.calculate_eic(local, MASS, FROM_RT, TO_RT, EIC_PPM, EIC_MZ)

    same = same_bits(eic_remote["x"], eic_local["x"]) and same_bits(
        eic_remote["y"], eic_local["y"]
    )

    print(f"file_bytes = {fixture.stat().st_size}")
    print(f"opening_bytes = {opening['bytes_sent']}")
    print(f"opening_requests = {opening['requests']}")
    print(f"query_bytes = {one_query['bytes_sent']}")
    print(f"query_requests = {one_query['requests']}")
    print(f"matches_local = {'yes' if same else 'no'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

import sys

MOST_REQUESTS_TO_OPEN = 12
MOST_REQUESTS_FOR_ONE_QUERY = 4
MOST_OF_THE_FILE_ONE_QUERY_MAY_TAKE = 0.1


def read(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            key, separator, value = line.strip().partition(" = ")
            if separator:
                values[key.strip()] = value.strip()
    return values


def main():
    measured = read(sys.argv[1])
    name = sys.argv[2]

    total = int(measured["file_bytes"])
    opening = int(measured["opening_bytes"])
    requests = int(measured["opening_requests"])
    query = int(measured["query_bytes"])
    requests_for_query = int(measured["query_requests"])
    matches = measured["matches_local"] == "yes"

    untouched = total - opening - query

    print(f"  file                 {total} bytes")
    print(f"  opening              {opening} bytes in {requests} requests")
    print(f"  one query            {query} bytes in {requests_for_query} requests")
    print(f"  never fetched        {untouched} bytes")

    failures = []

    if not matches:
        failures.append("the remote answer differs from the local one")

    if opening >= total:
        failures.append(
            f"opening fetched {opening} of {total} bytes: the whole file, so nothing is lazy"
        )

    if untouched <= 0:
        failures.append(
            f"every byte of the file was fetched: nothing was left, so nothing is lazy"
        )

    if query > total * MOST_OF_THE_FILE_ONE_QUERY_MAY_TAKE:
        failures.append(f"one query fetched {query} of {total} bytes: not lazy")

    if requests_for_query > MOST_REQUESTS_FOR_ONE_QUERY:
        failures.append(
            f"one query took {requests_for_query} requests: it is fetching range by range, not by plan"
        )

    if requests > MOST_REQUESTS_TO_OPEN:
        failures.append(
            f"opening took {requests} requests: the ranges are not being merged"
        )

    if failures:
        for line in failures:
            print(f"  FAIL  {line}")
        return 1

    share = 100 * untouched / total
    print(
        f"  ok  {name} is lazy: it never fetched {untouched} bytes, "
        f"{share:.0f}% of the file, and matches the local answer"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

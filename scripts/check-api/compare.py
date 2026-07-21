import struct
import sys


def read(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, separator, value = line.partition(" = ")
            if not separator:
                continue
            values[key.strip()] = value.strip()
    return values


def same(left, right):
    if left == right:
        return True
    try:
        return struct.pack(">d", float(left)) == struct.pack(">d", float(right))
    except ValueError:
        return False


def main():
    measured_path, reference_path, name, label = sys.argv[1:5]
    measured = read(measured_path)
    reference = read(reference_path)

    problems = []
    for key in sorted(set(measured) | set(reference)):
        if key not in measured:
            problems.append(f"    {key}: missing from {name}")
        elif key not in reference:
            problems.append(f"    {key}: {name} returned it, {label} does not have it")
        elif not same(measured[key], reference[key]):
            problems.append(
                f"    {key}: {name} got {measured[key]}, {label} says {reference[key]}"
            )

    if not problems:
        print(f"  ok  {name} matches {label}")
        return 0

    print(f"  FAIL  {name} differs from {label} in {len(problems)} places:")
    for line in problems[:20]:
        print(line)
    if len(problems) > 20:
        print(f"    ... and {len(problems) - 20} more")
    return 1


if __name__ == "__main__":
    sys.exit(main())

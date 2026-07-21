import os
import struct
import sys

import quantion

FROM_RT = 0.0
TO_RT = 5.0
EIC_PPM = 20.0
EIC_MZ = 0.005
GRID_START = 89.0
GRID_END = 91.0
GRID_STEP = 0.005
ROI_RANGE = 0.05
CORES = 1

FEATURE_MIN_INTENSITY = 500.0
FEATURE_MIN_PEAK_WIDTH_POINTS = 5

SELF_CHECK_VALUES = (1.0, 2.0, 3.0)

CHANNELS = [
    {
        "name": "mz_2",
        "mass": 90.05550,
        "min_intensity": 500.0,
        "min_peak_width_points": 5,
        "targets": [("met_1", 2.885), ("met_2", 2.552), ("met_3", 2.465)],
    },
    {
        "name": "mz_1",
        "mass": 89.04768,
        "min_intensity": 500.0,
        "min_peak_width_points": 5,
        "targets": [("met_1", 2.385), ("met_2", 2.18), ("met_3", 2.08)],
    },
]


def number(value):
    value = float(value)
    if value != value:
        return "nan"
    if value == float("inf"):
        return "inf"
    if value == float("-inf"):
        return "-inf"
    return repr(value)


def digest(values):
    hash_value = 0x811C9DC5
    for value in values:
        for byte in struct.pack("<d", float(value)):
            hash_value ^= byte
            hash_value = (hash_value * 0x01000193) & 0xFFFFFFFF
    return f"{hash_value:08x}"


def peak_options(channel):
    return {
        "min_intensity": channel["min_intensity"],
        "min_peak_width_points": channel["min_peak_width_points"],
        "auto_noise": True,
        "auto_baseline": True,
    }


def chrom_peak_options(channel):
    return {
        "min_intensity": channel["min_intensity"],
        "min_peak_width_points": channel["min_peak_width_points"],
        "auto_noise": False,
        "auto_baseline": False,
    }


def feature_peak_options():
    return {
        "min_intensity": FEATURE_MIN_INTENSITY,
        "min_peak_width_points": FEATURE_MIN_PEAK_WIDTH_POINTS,
        "auto_noise": True,
        "auto_baseline": True,
    }


def peak_lines(prefix, peak, out):
    shared_peak_lines(prefix, peak, out)
    out.append(f"{prefix}.n_points = {int(peak['n_points'])}")


def shared_peak_lines(prefix, peak, out):
    out.append(f"{prefix}.rt = {number(peak['rt'])}")
    out.append(f"{prefix}.from = {number(peak['from'])}")
    out.append(f"{prefix}.to = {number(peak['to'])}")
    out.append(f"{prefix}.integral = {number(peak['integral'])}")
    out.append(f"{prefix}.intensity = {number(peak['intensity'])}")


def chromatogram_index(channel):
    return sum(1 for other in CHANNELS if other["mass"] < channel["mass"])


def report_channel(file, channel, out):
    chrom_index = chromatogram_index(channel)
    name = channel["name"]
    eic = quantion.calculate_eic(
        file, channel["mass"], FROM_RT, TO_RT, EIC_PPM, EIC_MZ
    )
    x, y = list(eic["x"]), list(eic["y"])

    out.append(f"{name}.calculate_eic.points = {len(x)}")
    out.append(f"{name}.calculate_eic.x.fnv1a = {digest(x)}")
    out.append(f"{name}.calculate_eic.y.fnv1a = {digest(y)}")

    noise = quantion.find_noise_level(y)
    out.append(f"{name}.find_noise_level.intensity = {number(noise['intensity'])}")
    out.append(f"{name}.find_noise_level.width = {noise['width']}")

    found = quantion.find_peaks(x, y, peak_options(channel))
    out.append(f"{name}.find_peaks.count = {len(found)}")
    for index, peak in enumerate(found):
        peak_lines(f"{name}.find_peaks[{index}]", peak, out)

    target = channel["targets"][2][1]
    single = quantion.get_peak(x, y, target, ROI_RANGE, peak_options(channel))
    peak_lines(f"{name}.get_peak", single, out)

    targets = [
        {"id": id_, "rt": rt, "mz": channel["mass"], "range": ROI_RANGE}
        for id_, rt in channel["targets"]
    ]
    rows = quantion.get_peaks_from_eic(
        file, targets, FROM_RT, TO_RT, peak_options(channel), CORES
    )
    for row in rows:
        shared_peak_lines(f"{name}.get_peaks_from_eic.{row['id']}", row, out)

    items = [
        {"id": id_, "index": chrom_index, "rt": rt, "range": ROI_RANGE}
        for id_, rt in channel["targets"]
    ]
    chrom_rows = quantion.get_peaks_from_chrom(
        file, items, chrom_peak_options(channel), CORES
    )
    for (id_, _), row in zip(channel["targets"], chrom_rows):
        prefix = f"{name}.get_peaks_from_chrom.{id_}"
        out.append(f"{prefix}.rt = {number(row['rt'])}")
        out.append(f"{prefix}.from = {number(row['from'])}")
        out.append(f"{prefix}.to = {number(row['to'])}")
        out.append(f"{prefix}.intensity = {number(row['intensity'])}")
        out.append(f"{prefix}.area = {number(row['integral'])}")


def main():
    directory = sys.argv[1] if len(sys.argv) > 1 else "core/tests/fixtures/api"
    fixture = os.path.join(directory, "api.ion")

    file = quantion.parse_ion_path(fixture)
    out = []
    out.append(f"fnv1a.self_check = {digest(SELF_CHECK_VALUES)}")

    scans = quantion.get_scans(file, rt_range=(FROM_RT, TO_RT))
    out.append(f"parse_ion.scans = {len(scans)}")

    for channel in CHANNELS:
        report_channel(file, channel, out)

    grid = {"start": GRID_START, "end": GRID_END, "step_size": GRID_STEP}
    eic = {"ppm_tolerance": EIC_PPM, "mz_tolerance": EIC_MZ}

    found = quantion.find_features(
        file,
        FROM_RT,
        TO_RT,
        eic=eic,
        grid=grid,
        options=feature_peak_options(),
        cores=CORES,
    )
    out.append(f"find_features.count = {len(found)}")
    for index, feature in enumerate(found):
        name = f"find_features[{index}]"
        out.append(f"{name}.mz = {number(feature['mz'])}")
        out.append(f"{name}.rt = {number(feature['rt'])}")
        out.append(f"{name}.from = {number(feature['from'])}")
        out.append(f"{name}.to = {number(feature['to'])}")
        out.append(f"{name}.integral = {number(feature['integral'])}")
        out.append(f"{name}.intensity = {number(feature['intensity'])}")
        out.append(f"{name}.n_points = {int(feature['n_points'])}")

    consensus = quantion.get_features(
        directory,
        FROM_RT,
        TO_RT,
        eic=eic,
        grid=grid,
        options=feature_peak_options(),
        cores=CORES,
    )
    out.append(f"get_features.count = {len(consensus)}")
    for index, feature in enumerate(consensus):
        name = f"get_features[{index}]"
        out.append(f"{name}.mz = {number(feature['mz'])}")
        out.append(f"{name}.rt = {number(feature['rt'])}")
        out.append(f"{name}.from = {number(feature['from'])}")
        out.append(f"{name}.to = {number(feature['to'])}")
        out.append(f"{name}.integral = {number(feature['integral'])}")
        out.append(f"{name}.intensity = {number(feature['intensity'])}")

    print("\n".join(out))


if __name__ == "__main__":
    main()

import json
import os

here = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(here, "results.json")))

blue = "#2563eb"
orange = "#ea580c"
ink = "#111827"
grey = "#6b7280"
line = "#d1d5db"


def escape(text):
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def thousands(value):
    return f"{int(value):,}"


def bar_panel(x, width, title, axis_max, left_value, right_value, left_label, right_label, value_format):
    baseline = 330.0
    top = 110.0
    height = baseline - top
    bar_width = 62.0
    center = x + width / 2.0
    left_x = center - bar_width - 14.0
    right_x = center + 14.0

    def bar(bx, value, color):
        bar_height = 0.0 if axis_max <= 0 else value / axis_max * height
        by = baseline - bar_height
        label = value_format(value)
        return (
            f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" '
            f'rx="4" fill="{color}"/>'
            f'<text x="{bx + bar_width / 2:.1f}" y="{by - 10:.1f}" text-anchor="middle" '
            f'font-size="16" font-weight="700" fill="{ink}">{escape(label)}</text>'
        )

    return (
        f'<line x1="{x + 14:.1f}" y1="{baseline:.1f}" x2="{x + width - 14:.1f}" y2="{baseline:.1f}" '
        f'stroke="{line}" stroke-width="1.5"/>'
        + bar(left_x, left_value, blue)
        + bar(right_x, right_value, orange)
        + f'<text x="{left_x + bar_width / 2:.1f}" y="{baseline + 22:.1f}" text-anchor="middle" '
        f'font-size="13" fill="{grey}">{escape(left_label)}</text>'
        + f'<text x="{right_x + bar_width / 2:.1f}" y="{baseline + 22:.1f}" text-anchor="middle" '
        f'font-size="13" fill="{grey}">{escape(right_label)}</text>'
        + f'<text x="{center:.1f}" y="{baseline + 48:.1f}" text-anchor="middle" '
        f'font-size="15" font-weight="600" fill="{ink}">{escape(title)}</text>'
    )


def summary_figure():
    find = data["results"]["find_noise_level"]
    san = data["results"]["noise_san_plot"]
    width = 960
    height = 430
    panel_width = 300

    header = (
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>'
        f'<text x="{width / 2:.0f}" y="42" text-anchor="middle" font-size="24" font-weight="800" '
        f'fill="{ink}">find_noise_level vs noise_san_plot &#8212; TripleTOF, 970 truth features</text>'
        f'<text x="{width / 2:.0f}" y="70" text-anchor="middle" font-size="15" fill="{grey}">'
        f'same EICs, same settings, only the noise method differs &#183; higher recall is better, '
        f'fewer peaks is better</text>'
        f'<rect x="{width / 2 - 210:.0f}" y="84" width="14" height="14" rx="3" fill="{blue}"/>'
        f'<text x="{width / 2 - 190:.0f}" y="96" font-size="14" fill="{ink}">find_noise_level (current)</text>'
        f'<rect x="{width / 2 + 20:.0f}" y="84" width="14" height="14" rx="3" fill="{orange}"/>'
        f'<text x="{width / 2 + 40:.0f}" y="96" font-size="14" fill="{ink}">noise_san_plot</text>'
    )

    panels = (
        bar_panel(0, panel_width, "features found  (recall, /970)", 970,
                  find["features_found"], san["features_found"], "find", "san", lambda v: str(int(v)))
        + bar_panel(panel_width, panel_width, "sample detections  (/7760)", 8000,
                    find["sample_detections"], san["sample_detections"], "find", "san", thousands)
        + bar_panel(2 * panel_width, panel_width, "peaks returned  (junk, lower=better)", 110000,
                    find["peaks_returned"], san["peaks_returned"], "find", "san", thousands)
    )

    caption = (
        f'<text x="{width / 2:.0f}" y="{height - 12:.0f}" text-anchor="middle" font-size="15" '
        f'fill="{ink}">Recall is a tie (968 = 968). SanPlot returns '
        f'{data["derived"]["junk_ratio_sanplot_over_findnoise"]:.2f}x more peaks '
        f'(+{thousands(data["derived"]["extra_peaks_from_sanplot"])} junk).</text>'
    )

    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'font-family="Helvetica, Arial, sans-serif">'
        + header + panels + caption + "</svg>"
    )


def isolated_figure():
    cases = data["isolated_cases"]
    width = 720
    height = 400
    baseline = 300.0
    top = 110.0
    axis_max = 8
    chart_height = baseline - top
    group_width = 200.0
    bar_width = 54.0
    start_x = 60.0

    parts = [
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        f'<text x="{width / 2:.0f}" y="42" text-anchor="middle" font-size="22" font-weight="800" '
        f'fill="{ink}">Isolated-peak cases &#8212; samples detected (/8)</text>',
        f'<rect x="{width / 2 - 200:.0f}" y="60" width="14" height="14" rx="3" fill="{blue}"/>',
        f'<text x="{width / 2 - 180:.0f}" y="72" font-size="14" fill="{ink}">find_noise_level</text>',
        f'<rect x="{width / 2 + 20:.0f}" y="60" width="14" height="14" rx="3" fill="{orange}"/>',
        f'<text x="{width / 2 + 40:.0f}" y="72" font-size="14" fill="{ink}">noise_san_plot</text>',
        f'<line x1="40" y1="{baseline:.1f}" x2="{width - 40:.1f}" y2="{baseline:.1f}" '
        f'stroke="{line}" stroke-width="1.5"/>',
    ]

    for index, case in enumerate(cases):
        group_x = start_x + index * group_width
        center = group_x + group_width / 2.0
        left_x = center - bar_width - 8.0
        right_x = center + 8.0
        for bx, value, color in (
            (left_x, case["find_noise_level_samples"], blue),
            (right_x, case["san_plot_samples"], orange),
        ):
            bar_height = value / axis_max * chart_height
            by = baseline - bar_height
            parts.append(
                f'<rect x="{bx:.1f}" y="{by:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" '
                f'rx="4" fill="{color}"/>'
            )
            parts.append(
                f'<text x="{bx + bar_width / 2:.1f}" y="{by - 8:.1f}" text-anchor="middle" '
                f'font-size="15" font-weight="700" fill="{ink}">{int(value)}</text>'
            )
        parts.append(
            f'<text x="{center:.1f}" y="{baseline + 24:.1f}" text-anchor="middle" '
            f'font-size="14" font-weight="600" fill="{ink}">{escape(case["id"])}</text>'
        )
        parts.append(
            f'<text x="{center:.1f}" y="{baseline + 44:.1f}" text-anchor="middle" '
            f'font-size="12" fill="{grey}">mz {case["mz"]:.3f} &#183; rt {case["rt"]:.2f}</text>'
        )

    parts.append(
        f'<text x="{width / 2:.0f}" y="{height - 12:.0f}" text-anchor="middle" font-size="13" '
        f'fill="{grey}">TTOF0459 is past the 35 min window for both arms (not a noise effect).</text>'
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'font-family="Helvetica, Arial, sans-serif">' + "".join(parts) + "</svg>"
    )


open(os.path.join(here, "figure_summary.svg"), "w").write(summary_figure())
open(os.path.join(here, "figure_isolated.svg"), "w").write(isolated_figure())
print("wrote figure_summary.svg and figure_isolated.svg")

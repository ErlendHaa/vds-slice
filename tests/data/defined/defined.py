import segyio
import numpy as np
import sys

def create_defined(path):
    """ Create file with simple constant data.
    | xlines-ilines | 1                  | 2                  | 3                  |
    |---------------|--------------------|--------------------|--------------------|
    | 10            | 100, 101, 102, 103 | 108, 109, 110, 111 | 116, 117, 118, 119 |
    | 11            | 104, 105, 106, 107 | 112, 113, 114, 115 | 120, 121, 122, 123 |

    UTM coordinates for headers:
    | xlines-ilines | 1           | 2             | 3             |
    |---------------|-------------|---------------|---------------|
    | 10            | x=1, y=3    | x=2.1, y=3    | x=3.2, y=3    |
    | 11            | x=1, y=6.3  | x=2.1, y=6.3  | x=3.2, y=6.3  |
    """
    spec = segyio.spec()

    spec.sorting = 2
    spec.format = 1
    spec.samples = [0, 1, 2, 4]
    spec.ilines = [1, 2, 3]
    spec.xlines = [10, 11]

    # We use scaling constant of -10, meaning that values will be divided by 10
    il_step_x = int(1.1 * 10)
    il_step_y = int(0 * 10)
    xl_step_x = int(0 * 10)
    xl_step_y = int(3.3 * 10)
    ori_x = int(1 * 10)
    ori_y = int(3 * 10)

    with segyio.create(path, spec) as f:
        data = 100
        tr = 0
        for il in spec.ilines:
            for xl in spec.xlines:
                f.header[tr] = {
                    segyio.su.iline: il,
                    segyio.su.xline: xl,
                    segyio.su.cdpx:
                        (il - spec.ilines[0]) * il_step_x +
                        (xl - spec.xlines[0]) * xl_step_x +
                        ori_x,
                    segyio.su.cdpy:
                        (il - spec.ilines[0]) * il_step_y +
                        (xl - spec.xlines[0]) * xl_step_y +
                        ori_y,
                    segyio.su.scalco: -10,
                }
                data = data + len(spec.samples)
                f.trace[tr] = np.arange(start=data - len(spec.samples),
                                        stop=data, step=1, dtype=np.single)
                tr += 1

        f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING)


if __name__ == "__main__":
    path = sys.argv[1]
    create_defined(path)

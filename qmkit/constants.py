from pathlib import Path

HOME = Path(__file__).parent.resolve()

MEM_LIMIT = "2048 MB"

SAMPLE_MOLFILES = [str(p) for p in (HOME / "data").glob("*.mol")]


if __name__ == '__main__':
    print(SAMPLE_MOLFILES)

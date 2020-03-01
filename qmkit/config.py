from dataclasses import dataclass, asdict


CONFIG_LIST = {}


def register_config(cls):
    CONFIG_LIST[cls.__name__] = cls
    return cls


@register_config
@dataclass
class OptConfig:

    confgen: int = 1000

    rms: float = 1.0

    n_qm: int = 5


if __name__ == '__main__':
    config = CONFIG_LIST["OptConfig"]()
    print(asdict(config))
    print(config)

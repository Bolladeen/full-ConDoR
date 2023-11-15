def rgb_string_to_hex(rgb):
    """
    e.g. input: 'rgb(141,211,199)'
    """
    rgb = tuple(map(int, rgb[4:-1].split(',')))
    return '#%02x%02x%02x' % rgb

"""
Generate color map images and fill in the docs to point to them.

"""
from __future__ import print_function

import argparse
import os
import sys
from importlib import import_module

from jinja2 import Template
from palettable.palette import Palette

MODULES = {
    'palettable.colorbrewer.diverging': './colorbrewer/diverging',
    'palettable.colorbrewer.qualitative': './colorbrewer/qualitative',
    'palettable.colorbrewer.sequential': './colorbrewer/sequential',
    'palettable.tableau': './tableau',
    'palettable.wesanderson': './wesanderson',
    'palettable.cubehelix': './cubehelix'
}


def find_palettes(mod):
    """
    Find all Palette instances in mod.

    """
    return {
        k: v for k, v in vars(mod).items()
        if isinstance(v, Palette) and not k.endswith('_r')}


def gen_images(palettes, dir):
    """
    Create images for each palette in the palettes dict.
    For qualitative palettes only the discrete images is made.

    """
    img_dir = os.path.join(dir, 'img')
    os.makedirs(img_dir, exist_ok=True)

    discrete_fmt = '{}_discrete.png'.format
    continuous_fmt = '{}_continuous.png'.format

    img_size = (6, 0.5)

    for name, p in palettes.items():
        print('Making discrete image for palette {}'.format(name))
        p.save_discrete_image(
            os.path.join(img_dir, discrete_fmt(name)), size=img_size)

        if p.type != 'qualitative':
            print('Making continuous image for palette {}'.format(name))
            p.save_continuous_image(
                os.path.join(img_dir, continuous_fmt(name)), size=img_size)


def render_doc_page(dir, palette_names):
    """
    Render the documentation page in a given directory.

    """
    print('Rendering index in dir {}'.format(dir))

    with open(os.path.join(dir, 'index.md.tpl')) as f:
        tpl = Template(f.read())

    with open(os.path.join(dir, 'index.md'), 'w') as f:
        f.write(tpl.render(palettes=palette_names))


def mkdocs(mod, dir, images=False):
    palettes = find_palettes(mod)
    if images:
        gen_images(palettes, dir)
    render_doc_page(dir, sorted(palettes.keys()))


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description='Build Palettable documentation.')
    parser.add_argument(
        '-i', '--images', action='store_true', help='force rebuild images')
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    for mod, dir in MODULES.items():
        print('Running module {}'.format(mod))
        mkdocs(import_module(mod), dir, images=args.images)


if __name__ == '__main__':
    sys.exit(main())

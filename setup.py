from setuptools import setup, find_packages

setup(
    name='IntegrAlign',
    version='1.0',
    packages=find_packages(),
    install_requires=[
        'ttkthemes==2.1',
        'customtkinter',
        'napari[all]',
        'tk',
        'opencv-python',
        'SimpleITK',
        'scikit-image',
        'matplotlib',
        'imagecodecs',
    	'rasterio',
    	'shapely',
    	'geopandas',
    ],
    entry_points={
        'console_scripts': [
            'IntegrAlign = __main__:main'
        ]
    },
)

from setuptools import setup

setup(
    name='MSMuTect_3.2',
    version='3.2',
    packages=['src', 'src.Entry', 'src.GenomicUtils', 'src.IndelCalling', 'tests', 'tests.test_utils'],
    url='https://github.com/MaruvkaLab/MSMuTect_3.2',
    license='MIT',
    author='Avraham Kahan, Yossi Maruvka, Gaia Frant',
    author_email='avrahamkahan123@gmail.com',
    description='Tool to examine indels in tumors to determine microsatellite stability'
)

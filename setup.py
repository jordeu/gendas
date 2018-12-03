from distutils.core import setup

setup(
    name='gendas',
    version='0.1',
    packages=['gendas'],
    url='https://github.com/jordeu/gendas',
    license='Apache License 2.0',
    author='Jordi Deu-Pons',
    author_email='jordi@jordeu.net',
    description='Flexible and powerful genomic data manipulation library for Python',
    install_requires=['configobj', 'pathos', 'pytabix==0.0.2', 'bgdata', 'intervaltree', 'tqdm', 'click']
)

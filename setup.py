# coding=utf-8
import setuptools


setuptools.setup(name='zthash',
                 version='1.0.0',
                 description='Zémor-Tillich Hash Implementation',
                 long_description=open('README.md').read().strip(),
                 author='Itay Bookstein',
                 py_modules=['zthash'],
                 install_requires=['bitstring'],
                 license='MIT License',
                 zip_safe=False,
                 keywords='zthash')

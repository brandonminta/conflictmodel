from setuptools import setup, find_packages

setup(name='cfdadvection',
      description='Conflict simulation based on a tribute model',
      author='Brandon Minta',
      author_email='brandon.minta@yachaytech.edu.ec',
      packages=find_packages(),
      install_requires=['numpy', 'matplotlib', 'networkx', 'random'])
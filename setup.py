from setuptools import setup

setup(
   name='ppmp',
   version='1.0.0',
   description='Prediction of Perturbations of Modular Protein structures by triplet analysis',
   author='Vladimir Macko',
   author_email='vlamacko@gmail.com',
   packages=['ppmp'],
   install_requires=['matplotlib', 'numpy', 'pandas', 'seaborn', 'scipy', 'tqdm', 'biopython', 'statsmodels'],
)

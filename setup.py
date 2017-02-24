from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='kmbio',
      version='0.01.04',
      description='Cheap Personal bio-library',
      author='Keith Murray',
      author_email='kmurrayis@gmail.com',
      license='MIT',
      packages=['kmbio'],
      install_requires=[
          'numpy',
          'scikit-learn',
	  'ete3',
	  'biopython',
	  'scipy'
	  
      ],
      include_package_data=True,
      zip_safe=False)

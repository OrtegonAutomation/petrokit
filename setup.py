from setuptools import setup, find_packages

setup(
    name="petrokit",
    version="0.1.0",
    author="Camilo Andrés Ortegon Cuevas",
    author_email="camilo.ortegonc@outlook.com",
    description="Librería para análisis de transporte y producción en ingeniería de petróleos",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/coachito/petrokit",  # o tu repo cuando lo subas
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SROTuner",                            
    version="1.0",                            
    author="Yi Yao, Lin Li",
    author_email="yyao26@crimson.ua.edu",                        
    description="Python package to calculate, tune or explore the SRO in crystal structure.",
    long_description=long_description,          
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),        
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],                                         
    python_requires='>=3.10',                  
    py_modules=["SROTuner"],                        
    install_requires=[]                       
)
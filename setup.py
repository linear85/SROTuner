import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SROTuner",                            
    version="0.1",                            
    author="Yi Yao, Lin Li",                        
    description="A simple python package to tune the SRO of multi-component alloys",
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
    package_dir={'':'SROTuner/src'},            
    install_requires=[]                       
)
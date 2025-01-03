This is the supporting package for the projects in the papers listed below. 

Please feel free to reach out with errors, critiques, or questions.

Before the proper setup begins, it is assumed that you have conda installed and have run
conda init and restarted your shell. If you would like more information on this process
the link below is the official guide.

-

After this, open a blank shell and navigate to wherever you would like our code to be installed. Then use curl or git to pull our code. 

git clone linkname.git

Then you will need to use the afformentioned conda commands to build an enviroment file based on ours. 

conda env create --name envname --file=env.yml

Make sure to run this in the same directory as the yml file. This may take a while. 

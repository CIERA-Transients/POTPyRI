procedure disp(inlist,reference)

# Script to dispersion correct a list a files with a given reference image
# all it does is elimate the annoyance of having to hedit the
# file before doing a dispcor
# 05/2001  RC

string	inlist		{prompt="Image(s) for dispersion correction"}
string	reference	{prompt="Reference image"}
string	prefix		{"d",prompt="Prefix for output images"}

struct	*inimglist

begin

	string 	infile, img

	infile =  mktemp("tmp$lick")
	sections (inlist,option="fullname",>infile)
	inimglist = infile
	while (fscan(inimglist,img) != EOF) {
	      hedit(img,"REFSPEC1",reference,del+,add+,ver-,show+)
	      dispcor(img,prefix // img)
	}

	inimglist = ""

end

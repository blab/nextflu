var	vaccineChoice = {
    "B/Beijing/184/1993": "1998-11-01",
    "B/Sichuan/379/1999": "2001-09-25",
    "B/Shanghai/361/2002": "2004-09-25",
    "B/Florida/4/2006": "2008-09-25",
    "B/Wisconsin/1/2010": "2012-02-25",
    "B/Massachusetts/2/2012": "2013-02-25",
    "B/Phuket/3073/2013": "2014-09-25"
};
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]};
var default_gene = 'HA1';

var structure = "4M40.pdb";

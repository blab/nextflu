var	vaccineChoice = {};
vaccineChoice['B/Beijing/184/93'] = "1998-11-01";
vaccineChoice['B/Sichuan/379/99'] = "2001-09-25";
vaccineChoice['B/Shanghai/361/2002'] = "2004-09-25";
vaccineChoice['B/Florida/4/2006'] = "2008-09-25";
vaccineChoice['B/Wisconsin/01/2010'] = "2012-02-25";
vaccineChoice['B/Massachusetts/02/2012'] = "2013-02-25";
vaccineChoice['B/PHUKET/3073/2013'] = "2014-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]}
var structure = "4M40.pdb"

var	vaccineChoice = {
    "B/Shangdong/7/1997": "1999-09-25",
    "B/HongKong/330/2001": "2002-09-25",
    "B/Malaysia/2506/2004": "2006-09-25",
    "B/Brisbane/60/2008": "2009-09-25",
    "B/Colorado/6/2017": "2018-02-22",
    "B/Washington/2/2019": "2019-09-27"
};
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]};
var default_gene = 'HA1';

var structure = "4FQM.pdb";

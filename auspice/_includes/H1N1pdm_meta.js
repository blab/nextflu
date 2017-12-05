var	vaccineChoice = {};
vaccineChoice['A/California/7/2009'] = "2009-09-25";
vaccineChoice['A/Michigan/45/2015'] = "2016-09-29";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [1,20,52]],
                         'HA1':[[1,1,1], [52,460,52+981]],
                         'HA2':[[1.2,1.2,1.2], [52+981,1200,1701]]};
var default_gene = 'HA1';

var structure = "4LXV.pdb"

var	vaccineChoice = {};
vaccineChoice['A/California/07/2009'] = "2009-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [1,20,52]],
                         'HA1':[[1,1,1], [52,460,52+981]],
                         'HA2':[[1.2,1.2,1.2], [52+981,1200,1701]]}
var structure = "4LXV.pdb"

var reference_viruses = {
    "A/England/195/2009":'2009-04-28',
    "A/Auckland/3/2009":'2009-04-25',
    "A/HongKong/2212/2010":'2010-07-16',
    "A/England/106/2010":'2010-05-01',
    "A/Lviv/N6/2009":'2009-10-27',
    "A/HongKong/2200/2010":'2010-07-14',
    "A/Bayern/69/2009":'2009-05-01',
    "A/Astrakhan/1/2011":'2011-02-28',
    "A/SouthAfrica/3626/2013":'2013-06-06',
    "A/England/87/2010":'2010-01-01',
    "A/Narita/1/2009":'2009-05-08',
    "A/HongKong/3934/2011":'2011-03-29',
    "A/HongKong/5659/2012":'2012-05-21'};

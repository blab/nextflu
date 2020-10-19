var	vaccineChoice = {
    "A/California/1/2009": "2009-09-25",
    "A/Michigan/45/2015": "2016-09-29",
    "A/Brisbane/2/2018": "2019-02-20",
    "A/Guangdong-Maonan/SWL1536/2019-egg": "2020-02-28",
    "A/Hawaii/70/2019": "2020-02-28",
    "A/Victoria/2570/2019-egg": "2020-09-25",
    "A/Wisconsin/588/2019": "2020-09-25"
};

var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [1,20,52]],
                         'HA1':[[1,1,1], [52,460,52+981]],
                         'HA2':[[1.2,1.2,1.2], [52+981,1200,1701]]};
var default_gene = 'HA1';

var structure = "4LXV.pdb"

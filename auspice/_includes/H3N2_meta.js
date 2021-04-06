var	vaccineChoice = {
    "A/Sydney/5/1997": "1997-09-25",
    "A/Moscow/10/1999": "1999-09-25",
    "A/Fujian/411/2002": "2003-09-25",
    "A/California/7/2004": "2005-02-21",
    "A/Wisconsin/67/2005": "2006-02-21",
    "A/Brisbane/10/2007": "2007-09-25",
    "A/Perth/16/2009": "2009-09-25",
    "A/Victoria/361/2011": "2012-02-21",
    "A/Texas/50/2012": "2013-09-25",
    "A/Switzerland/9715293/2013": "2014-09-25",
    "A/HongKong/4801/2014": "2015-09-24",
    "A/Singapore/Infimh-16-0019/2016": "2017-09-28",
    "A/Switzerland/8060/2017": "2018-09-27",
    "A/Kansas/14/2017": "2019-03-21",
    "A/SouthAustralia/34/2019": "2019-09-27",
    "A/HongKong/2671/2019-egg": "2020-02-28",
    "A/HongKong/45/2019": "2020-02-28",
    "A/Cambodia/e0826360/2020": "2021-02-26"
};
var vaccineStrains = Object.keys(vaccineChoice);
var branch_labels=false;
var restrictTo = {"region": "all"};

var genome_annotation = {'SP': [[1.2,1.2,1.2], [1,20,49]],
                         'HA1': [[1,1,1], [49,460,49+987]],
						 'HA2': [[1.2,1.2,1.2], [49+987,1200,1701]]};
var default_gene = 'HA1';

var structure = "5HMG.pdb"

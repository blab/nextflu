var	vaccineChoice = {};
vaccineChoice['B/Shangdong/7/97'] = "1999-09-25";
vaccineChoice['B/HongKong/330/2001'] = "2002-09-25";
vaccineChoice['B/Malaysia/2506/2004'] = "2006-09-25";
vaccineChoice['B/Brisbane/60/2008'] = "2009-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'HA1':[[1,1,1], [1,160,346]],
                         'HA2':[[1.2,1.2,1.2], [346,450,600]]}

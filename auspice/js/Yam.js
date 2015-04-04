var file_prefix = 'Yam_';
var dfreq_dn = 2;

var	vaccineChoice = {};
vaccineChoice['B/Beijing/184/93'] = "1998-11-01";
vaccineChoice['B/Sichuan/379/99'] = "2001-09-25";
vaccineChoice['B/Shanghai/361/2002'] = "2004-09-25";
vaccineChoice['B/Florida/4/2006'] = "2008-09-25";
vaccineChoice['B/Wisconsin/01/2010'] = "2012-02-25";
vaccineChoice['B/Massachusetts/02/2012'] = "2013-02-25";
vaccineChoice['B/PHUKET/3073/2013'] = "2014-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var epiColorDomain = [0,1,2,3, 4,5,6,7,8,9,10];
var nonEpiColorDomain = [0,1,2,3,4,5,6,7,8,9,10,11];
var rbsColorDomain = [0,0.5, 1,1.5, 2];

var time_ticks = [2010, 2011, 2012, 2013, 2014, 2015];
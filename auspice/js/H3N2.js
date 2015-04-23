var file_prefix = 'H3N2_HI_';
var dfreq_dn = 2;

var	vaccineChoice = {};
vaccineChoice['A/Fujian/411/2002'] = "2003-09-25";
vaccineChoice['A/California/7/2004'] = "2005-02-21";
vaccineChoice['A/Wisconsin/67/2005'] = "2006-02-21";
vaccineChoice['A/Brisbane/10/2007'] = "2007-09-25";
vaccineChoice['A/Perth/16/2009'] = "2009-09-25";
vaccineChoice['A/Victoria/361/2011'] = "2012-02-21";
vaccineChoice['A/Texas/50/2012'] = "2013-09-25";
vaccineChoice['A/Switzerland/9715293/2013'] = "2014-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var epiColorDomain = [4,5,6,7,8,9,10,11,12,13];
var nonEpiColorDomain = [2,3,4,5,6,7,8,9,10,11];
var rbsColorDomain = [0,1,2,3,4];
var dfreqColorDomain = [-0.18, -0.14, -0.10, -0.06, -0.02, 0.02, 0.06, 0.10, 0.14, 0.18];
var genericCD = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
var HIColorDomain = genericCD.map(function (d){return 15*d;});

var time_ticks = [2012.5, 2013, 2013.5, 2014, 2014.5, 2015];
var	time_window = 1.0;  // layer of one year that is considered current or active
var LBItau = 0.0008;
var LBItime_window = 0.5;


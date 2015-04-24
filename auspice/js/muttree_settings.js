var file_prefix = '';
var path = '';
var dfreq_dn = 2;
var branch_labels = true;
var tip_labels = true;

var	vaccineChoice = {};
var vaccineStrains = Object.keys(vaccineChoice);

var epiColorDomain = [4,5,6,7,8,9,10,11,12,13];
var nonEpiColorDomain = [2,3,4,5,6,7,8,9,10,11];
var dateDomain = [2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013];
var rbsColorDomain = [0,1,2,3,4];
var dfreqColorDomain = [-0.18, -0.14, -0.10, -0.06, -0.02, 0.02, 0.06, 0.10, 0.14, 0.18];

var time_ticks = [2012.5, 2013, 2013.5, 2014, 2014.5, 2015];
var	time_window = 1.0;  // layer of one year that is considered current or active
var LBItau = 0.0008;
var LBItime_window = 0.5;


var	vaccineChoice = {};
vaccineChoice['B/Beijing/184/1993'] = "1998-11-01";
vaccineChoice['B/Sichuan/379/1999'] = "2001-09-25";
vaccineChoice['B/Shanghai/361/2002'] = "2004-09-25";
vaccineChoice['B/Florida/4/2006'] = "2008-09-25";
vaccineChoice['B/Wisconsin/1/2010'] = "2012-02-25";
vaccineChoice['B/Massachusetts/2/2012'] = "2013-02-25";
vaccineChoice['B/Phuket/3073/2013'] = "2014-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]};
var default_gene = 'HA1';

var structure = "4M40.pdb";
var reference_viruses = {
"B/Florida/7/2004":'2004-04-01',
"B/Brisbane/3/2007":'2007-03-09',
"B/Serbia/1894/2011":'2011-03-08',
"B/Valladolid/18/2008":'2008-02-18',
"B/Novosibirsk/1/2012":'2012-02-14',
"B/Algeria/G486/2010":'2010-06-06',
"B/Stockholm/12/2011":'2011-02-28',
"B/Niedersachsen/1/2010":'2010-10-18',
"B/Utah/9/2014":'2014-05-29',
"B/Massachusetts/2/2012":'2012-03-13',
"B/HongKong/3577/2012":'2012-06-13',
"B/Phuket/3073/2013":'2013-11-21',
"B/Florida/4/2006":'2006-11-01',
"B/Estonia/55669/2011":'2011-03-14',
"B/England/145/2008":'2008-05-01'
};

var countries = [];
var divisions = [];
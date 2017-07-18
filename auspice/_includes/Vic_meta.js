var	vaccineChoice = {};
vaccineChoice['B/Shangdong/7/97'] = "1999-09-25";
vaccineChoice['B/HongKong/330/2001'] = "2002-09-25";
vaccineChoice['B/Malaysia/2506/2004'] = "2006-09-25";
vaccineChoice['B/Brisbane/60/2008'] = "2009-09-25";
var vaccineStrains = Object.keys(vaccineChoice);

var restrictTo = {"region":"all"};

var genome_annotation = {'SP':[[1.2,1.2,1.2], [12,20,57]],
                         'HA1':[[1,1,1], [57,460,57+1038]],
                         'HA2':[[1.2,1.2,1.2], [57+1038,1200,1769]]};
var default_gene = 'HA1';

var structure = "4FQM.pdb";

var reference_viruses = {
"B/Odessa/3886/2010":'2010-03-19',
"B/Victoria/304/2006":'2006-06-06',
"B/HongKong/22/2001":'2001-01-18',
"B/Brisbane/60/2008":'2008-08-04',
"B/Madagascar/7002/2009":'2009-08-31',
"B/Beijing/243/1997":'1997-01-01',
"B/Akita/5/2001":'2001-08-01',
"B/HongKong/514/2009":'2009-10-11',
"B/Cambodia/30/2011":'2011-01-04',
"B/Aichi/3/1998":'1998-11-01',
"B/HongKong/45/2005":'2005-02-07',
"B/Malaysia/2506/2004":'2004-07-01',
"B/Hawaii/33/2004":'2004-03-01',
"B/Brazil/2937/2008":'2008-07-28',
"B/Texas/26/2008":'2008-11-24',
"B/HongKong/330/2001":'2001-07-01',
"B/HongKong/335/2001":'2001-11-01',
"B/Paris/1762/2009":'2009-02-09',
"B/SouthAustralia/81/2012":'2012-11-28',
"B/Formosa/V2367/2012":'2012-08-06',
"B/HongKong/1434/2002":'2002-01-18',
"B/Brisbane/32/2002":'2002-03-01',
"B/Akita/27/2001":'2001-06-01',
"B/Shandong/7/1997":'1997-12-01',
"B/Brisbane/33/2008":'2008-07-13',
"B/NewCaledonia/5/2006":'2006-04-24'
};

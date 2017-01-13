var frequencies, pivots, entropy;
var gene = 'nuc';
var mutType = 'aa';
var plot_frequencies = true;

/**
 * for each node, calculate the derivative of the frequency tranjectory. if none exists, copy parent
**/
function calcDfreq(node, freq_ii){
	if (typeof node.children != "undefined") {
		for (var i1=0; i1<node.children.length; i1++) {
			var label_str = "global_clade:"+node.children[i1].clade;
			if (typeof frequencies !== "undefined" && frequencies[label_str] != "undefined"){
				var tmp_freq = get_frequencies("global", "clade:"+node.children[i1].clade)
				node.children[i1].dfreq = (tmp_freq[freq_ii] + 0.01)/(tmp_freq[freq_ii-dfreq_dn] + 0.01);
			} else {
				node.children[i1].dfreq = node.dfreq;
			}
			calcDfreq(node.children[i1], freq_ii);
		}
	}
};

/**
loops over all genotypes from a certain region and sums the frequency contributions
of the genotype matches at the specified positions
**/
function get_frequencies(region, gt){
	var freq = [];
	for (var pi=0; pi<pivots.length; pi++){freq[freq.length]=0;}
	console.log("searching for "+region+' ' + gt);
	var label_str = region+'_'+gt
	if (frequencies[label_str]!=undefined) {
		console.log(gt+" found");
		for (var pi=0; pi<freq.length; pi++){
			freq[pi]+=frequencies[label_str][pi];
		}
	}else{
		console.log("not found "+gt);
	}
	return freq.map(function (d) {return Math.round(d*100)/100;});
};

var freqDataString = "";
function make_gt_chart(gt){
	var tmp_data = [];
	var tmp_trace = ['x'];
	var tmp_colors = {};
	tmp_data.push(tmp_trace.concat(pivots));
	gt.forEach(function (d, i) {
		var region = d[0];
		var genotype = d[1];
		var freq = get_frequencies(region, genotype);
		console.log(region+' '+genotype);
		if (d3.max(freq)>0) {
			var tmp_trace = genotype.toString().replace(/,/g, ', ');
			if (region != "global") {
				tmp_trace = region + ':\t' + tmp_trace;
			}
			tmp_data.push([tmp_trace].concat(freq));
			tmp_colors[tmp_trace] = genotypeColors[i];
		}
	});
	console.log(tmp_colors);
	gt_chart.load({
       	columns: tmp_data,
       	unload: true
	});
	gt_chart.data.colors(tmp_colors);
	// construct a tab separated string the frequency data
	freqDataString="";
	for (var ii=0; ii<tmp_data[0].length; ii+=1){
		for (var jj=0; jj<tmp_data.length; jj+=1){
			freqDataString += "" + tmp_data[jj][ii] + ((jj<tmp_data.length-1)?"\t":"\n");
		}
	}
}

function addClade(d) {
	if (typeof gt_chart != "undefined"){
		console.log(d);
		var plot_data = [['x'].concat(pivots)];
		var reg = "global";
		var label_str = 'clade:'+d.target.clade;
		if (typeof frequencies[reg+"_"+label_str] !="undefined" ){

			plot_data[plot_data.length] = [reg].concat(get_frequencies(reg,label_str));
		}
		if (plot_data.length > 1) {
			if (plot_data[1][0] == reg) {
				plot_data[1][0] = "clade";
			}
		}
		gt_chart.load({
	       	columns: plot_data
		});
	}
}

function removeClade() {
	if (typeof gt_chart != "undefined"){
		gt_chart.unload({
	       	ids: ["clade"]
		});
	}
}

width = parseInt(d3.select(".freqplot-container").style("width"), 10);
height = 250;
var position = "inset";

var gt_chart = c3.generate({
	bindto: '#gtchart',
	size: {width: width-10, height: height},
	onresize: function() {
		width = parseInt(d3.select(".freqplot-container").style("width"), 10);
		height = 250;
		gt_chart.resize({height: height, width: width});
	},
	legend: {
		position: position,
		inset: {
    		anchor: 'top-right',
    		x: 10,
    		y: -15,
    		step: 1
    	}
	},
  	color: {
        pattern: genotypeColors
    },
	axis: {
		y: {
			label: {
				text: 'frequency',
				position: 'outer-middle'
			},
			tick: {
				values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
				outer: false
			},
            min: 0,
			max: 1
		},
		x: {
			label: {
				text: 'time',
				position: 'outer-center'
			},
			tick: {
				values: time_ticks,
				outer: false
			}
		}
	},
	data: {
		x: 'x',
		columns: [],
	}
});

d3.json(path + file_prefix + "frequencies.json", function(error, json){
	console.log(error, path+file_prefix+"frequencies.json");
	frequencies = json;
	console.log(frequencies);
	pivots= frequencies["pivots"].map(function (d) {return Math.round(parseFloat(d)*100)/100;});
	var ticks = [Math.round(pivots[0])];
	var tick_step = Math.round((pivots[pivots.length-1]-pivots[0])/6*10)/10;
	while (ticks[ticks.length-1]<pivots[pivots.length-1]){
		ticks.push(Math.round((ticks[ticks.length-1]+tick_step)*10)/10);
	}

	d3.select("#plotfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);
			make_gt_chart(gt);
		});
	d3.select("#downloadfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);
			make_gt_chart(gt);
			var blob = new Blob([freqDataString], {type: "text/plain;charset=utf-8"});
			saveAs(blob,'frequencies.tsv');
		});
	make_gt_chart(parse_gt_string(document.getElementById("gtspec").value));
});

/**
 * for each node, calculate the derivative of the frequency tranjectory. if none exists, copy parent
**/
function calcDfreq(node, freq_ii){
	if (typeof node.children != "undefined") {
		for (var i1=0; i1<node.children.length; i1++) {
			if (typeof node.children[i1].freq != "undefined") {
				if (node.children[i1].freq["global"] != "undefined"){
					var tmp_freq = node.children[i1].freq["global"]
					node.children[i1].dfreq = 0.5*(tmp_freq[freq_ii] - tmp_freq[freq_ii-dfreq_dn])/(tmp_freq[freq_ii] + tmp_freq[freq_ii-dfreq_dn] + 0.1);
				} else {
					node.children[i1].dfreq = node.dfreq;
				}
			}
			calcDfreq(node.children[i1], freq_ii);
		}
	}
};


/**
 * for each node, calculate the number of tips in the currently selected time window. 
**/
function calcTipCounts(node){
	node.tipCount = 0;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
			calcTipCounts(node.children[i]);
			node.tipCount += node.children[i].tipCount;
		}
	}
	else if (node.current){ 
		node.tipCount = 1;
	}
};

width = parseInt(d3.select(".freqplot-container").style("width"), 10);
var position = "right";
if (width < 600) {
	position = "bottom";
}

var gt_chart = c3.generate({
	bindto: '#gtchart',
	size: {width: width, height: 350},
	legend: {position: position},
  	color: {
        pattern: ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]
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
		columns: []
	}
});



function adjust_freq_by_date() {
	calcTipCounts(rootNode);
	var tipCount = rootNode.tipCount;		
	console.log("Total tipcount: " + tipCount);	
	nodes.forEach(function (d) {
		d.frequency = (d.tipCount)/tipCount;
	});
}	

function contains(arr, obj) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == obj) return true;
    }
}

d3.json("/data/" + file_prefix + "frequencies.json", function(error, json){
	console.log(error);
	var pivots= json["mutations"]["global"]["pivots"].map(function (d) {return Math.round(parseFloat(d)*100)/100;});
	var ticks = [Math.round(pivots[0])];
	var step = Math.round((pivots[pivots.length-1]-pivots[0])/6*10)/10;
	while (ticks[ticks.length-1]<pivots[pivots.length-1]){
		ticks.push(Math.round((ticks[ticks.length-1]+step)*10)/10);
	}
	//gt_chart.axis.x.values = ticks;
	/**
		parses a genotype string into region and positions
	**/
	function parse_gt_string(gt){
		separate_plots = gt.split(',');
		mutations = separate_plots.map(
			function (d) {	var tmp = d.split(/[\s//]/); //FIXME: make more inclusive
							var region;
							var positions = [];
							for (var i=0; i<tmp.length; i++){
								if (contains(["EU","NA","AS","OC"], tmp[i])){
									region = tmp[i];
								}else{
									if (tmp[i].length>0) positions.push(tmp[i]);
								}
							}
							if (typeof region == "undefined") region="global";
							// sort of this is a multi mutation genotype
							if (positions.length>1){
								positions.sort(function (a,b){
									return parseInt(a.substring(0,a.length-1)) - parseInt(b.substring(0,b.length-1));
								});
							}
							return [region, positions.join('/')];});
		return mutations;
	};

	/**
	loops over all genotypes from a certain region and sums the frequency contributions
	of the genotype matches at the specified positions
	**/
	function get_frequencies(region, gt){
		var freq = [];
		for (var pi=0; pi<pivots.length; pi++){freq[freq.length]=0;}
		if (json["clades"][region][gt.toLowerCase()]!=undefined) {
			console.log(gt+" found as clade");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=json["clades"][region][gt.toLowerCase()][pi];
			}
		}
		else if ((typeof json["genotypes"] !="undefined") && (json["genotypes"][region][gt]!=undefined)) {
			console.log(gt+" found as genotype");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=json["genotypes"][region][gt][pi];
			}
		}else if (json["mutations"][region][gt]!=undefined) {
			console.log(gt+" found as mutation");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=json["mutations"][region][gt][pi];
			}
		}
		return freq.map(function (d) {return Math.round(d*100)/100;});
	};

	function make_gt_chart(gt){
		var tmp_data = [];
		var tmp_trace = ['x'];
		tmp_data.push(tmp_trace.concat(pivots));
		gt.forEach(function (d) {
			var region = d[0];
			var genotype = d[1];
			var freq = get_frequencies(region, genotype);
			var tmp_trace = genotype.toString().replace(/,/g, ', ');
			if (region != "global") {
				tmp_trace = region + ':\t' + tmp_trace;
			}
			tmp_data.push([tmp_trace].concat(freq));
		});
		gt_chart.load({
	       	columns: tmp_data,
	       	unload: true
		});
	}

	d3.select("#plotfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);			
			make_gt_chart(gt);
		});
	make_gt_chart(parse_gt_string(document.getElementById("gtspec").value));
});

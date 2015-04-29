var frequencies, pivots;

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

function contains(arr, obj) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == obj) return true;
    }
}

d3.json(path + file_prefix + "frequencies.json", function(error, json){
	console.log(error);
	frequencies = json;
	pivots= frequencies["mutations"]["global"]["pivots"].map(function (d) {return Math.round(parseFloat(d)*100)/100;});
	var ticks = [Math.round(pivots[0])];
	var tick_step = Math.round((pivots[pivots.length-1]-pivots[0])/6*10)/10;
	while (ticks[ticks.length-1]<pivots[pivots.length-1]){
		ticks.push(Math.round((ticks[ticks.length-1]+tick_step)*10)/10);
	}
	//gt_chart.axis.x.values = ticks;
	/**
		parses a genotype string into region and positions
	**/

	var ent = [['x'],['entropy']];
	for (var ii=0;ii<frequencies["entropy"].length;ii+=1){
		if (Math.round(10000*frequencies["entropy"][ii][1])/10000>0){
			ent[1].push(Math.round(10000*frequencies["entropy"][ii][1])/10000);
			ent[0].push(ii+1);
		}
	}

	
	var entropy_chart = c3.generate({
		bindto: '#entropy',
		size: {width: width, height: 350},
		legend: {position: position},
		axis: {
			y: {
				label: {
					text: 'entropy',
					position: 'outer-middle'	
				}
			},
			x: {
				label: {
					text: 'position',
					position: 'outer-right'	
				},
				tick: {
					outer: false,
					values: [100,200,300,400,500]				
				}
			}
		},			
		data: {
			x: 'x',
			columns: ent,
			type: "bar",
			onclick: function (d,i) { 
						console.log(d);
		            	if (frequencies["entropy"][d.x-1][2].length>1){
		            		var tmp = [];
		            		for (var ii=0;ii<frequencies["entropy"][d.x-1][2].length;ii+=1){
								tmp.push(["global",d.x+frequencies["entropy"][d.x-1][2][ii]]);
		            		}
		            		console.log(tmp);
		            		make_gt_chart(tmp);
		            		colorBy = "genotype";
		            		colorByGenotypePosition([d.x-1]);
		            	}
		            }
		},
	    tooltip: {
	        format: {
	            title: function (d) { 
	            	return 'Position ' + d + frequencies["entropy"][d-1][2].join(","); },
	            value: function (value, ratio, id) {
	                return "Entropy: "+value;
	            }
	        }
		},
	});

	d3.select("#plotfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);			
			make_gt_chart(gt);
		});
	make_gt_chart(parse_gt_string(document.getElementById("gtspec").value));
});

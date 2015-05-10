var regions = ["Africa", "SouthAmerica", "WestAsia", "Oceania", "Europe", "JapanKorea", "NorthAmerica", "SoutheastAsia", "SouthAsia", "China"]
var restrictTo = "all";

var cladeToSeq = {}

var globalDate = new Date();

var nodes, tips, rootNode, links, vaccines;

var nDisplayTips, displayRoot;

function treePlotHeight(width) {
	return 400 + 0.35*width;
}
var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
var treeWidth = containerWidth;
var treeHeight = treePlotHeight(treeWidth);
var tree = d3.layout.tree()
	.size([treeHeight, treeWidth]);


var treeplot = d3.select("#treeplot")
	.attr("width", treeWidth)
	.attr("height", treeHeight);

var legend = d3.select("#legend")
	.attr("width", 280)
	.attr("height", 100);

var colorBy = document.getElementById("coloring").value;
var colorScale;

var time_step;


d3.json(path + file_prefix + "meta.json", function(error, json) {
    if (error) return console.warn(error);
    d3.select("#updated").text(json['updated']);
    commit_id = json['commit'];
    short_id = commit_id.substring(0, 6);   
    d3.select("#commit")
        .append("a")
        .attr("href", "http://github.com/blab/nextflu/commit/" + commit_id)
        .text(short_id);


    if (typeof json['virus_stats_before_subsampling'] != "undefined"){
        width = parseInt(d3.select(".virus_stats-container").style("width"), 10);
        height = 250;
        var tmp_trace = ['x'];
        var available_viruses_count = [];
        available_viruses_count.push(tmp_trace.concat(json['dates']));
        for (var i=0; i<json['regions'].length; i++){
            reg = json['regions'][i];
            tmp_trace = [reg];
            available_viruses_count.push(tmp_trace.concat(json['virus_stats_before_subsampling'][reg].map(function(d) {return Math.sqrt(Math.max(d+0.0));})));
        }
        var presub_virus_chart = c3.generate({
            bindto: '#available_viruses',
            size: {width: width-10, height: height},
            onresize: function() {
                width = parseInt(d3.select(".available_viruses").style("width"), 10);
                height = 250;
                presub_virus_chart.resize({height: height, width: width});
            },      
            color: {pattern: regionColors},
            legend: {show: true},
            axis: {
                y: {
                    label: {
                        text: '# of available sequences',
                        position: 'outer-middle'
                    },
                    tick: {
                        //format: function (d) { return Math.pow(10,d).toFixed(0); }
                        format: function (d) { return (d*d).toFixed(0); },
                        values: [0,1,5,10,20]
                    }           
                },
                x: {
                    label: {
                        text: 'time',
                        position: 'outer-center',
                    },
                    tick: {
                        outer: false,
                        values: time_ticks
                    }
                }
            },          
            data: {
                x: 'x',
                columns: available_viruses_count,
            },
            tooltip: {
                format: {
                    value: function (d) { return (d*d).toFixed(0);},
                    title: function (d) {return "Date: "+d.toFixed(2);}
                }
            }
        });
    }
    if (typeof json['virus_stats'] != "undefined"){
        var tmp_trace = ['x'];
        var sampled_virus_count = [];
        sampled_virus_count.push(tmp_trace.concat(json['dates']));
        var all_regions = ['all'];
        for (var i=0; i<sampled_virus_count[0].length-1; i++){all_regions.push(0);}
        for (var i=0; i<json['regions'].length; i++){
            reg = json['regions'][i];
            tmp_trace = [reg];
            sampled_virus_count.push(tmp_trace.concat(json['virus_stats'][reg]));
            for (var j=0; j<json['virus_stats'][reg].length; j++){all_regions[j+1]+=json['virus_stats'][reg][j];}
        }
        sampled_virus_count.push(all_regions);
        var virus_stats_chart = c3.generate({
            bindto: '#sampled_viruses',
            size: {width: width-10, height: height},
            onresize: function() {
                width = parseInt(d3.select(".sampled_viruses").style("width"), 10);
                height = 250;
                virus_stats_chart.resize({height: height, width: width});
            },      
            color: {pattern: regionColors},
            axis: {
                y: {
                    label: {
                        text: '# of sampled sequences',
                        position: 'outer-middle'
                    },
                    tick: {
                        values: [0,10,20],
                    }           
                },
                x: {
                    label: {
                        text: 'time',
                        position: 'outer-center',
                    },
                    tick: {
                        outer: false,
                        values: time_ticks
                    }
                }
            },          
            data: {
                x: 'x',
                columns: sampled_virus_count,
            },
            tooltip: {
                format: {
                    title: function (d) {return "Date: "+d.toFixed(2);}
                }
            }
        });
    }
});

var legendRectSize = 15;
var legendSpacing = 4;
function makeLegend(){
	
	d3.select("#legend-title").text(function(d){
		if (colorBy == "ep") {
			return "Epitope mutations";
		}
		if (colorBy == "ne") {
			return "Non-epitope mutations";
		}
		if (colorBy == "rb") {
			return "Receptor binding mutations";
		}
		if (colorBy == "lbi") {
			return "Local branching index";
		}
		if (colorBy == "region") {
			return "Region";
		}
		if (colorBy == "genotype") {
			return "Genotype";
		}
		if (colorBy == "date") {
			return "Date";
		}		
		if (colorBy == "dfreq") {
			var tmp_nmonth = Math.round(12*dfreq_dn*time_step);
			var tmp_text = "Freq. change ("+tmp_nmonth+" month";
			if (tmp_nmonth>1){
				tmp_text+='s';
			}
			return tmp_text+')';
		}
    });

    if (colorBy == "region"){
        make_map();
    }else{
      var tmp_leg = legend.selectAll(".legend")
      .data(colorScale.domain())
      .enter().append('g')
      .attr('class', 'legend')
      .attr('transform', function(d, i) {
       var stack = 5;
       var height = legendRectSize + legendSpacing;
       var fromRight = Math.floor(i / stack);
       var fromTop = i % stack;
       var horz = fromRight * 145 + 5;				
       var vert = fromTop * height + 5;
       return 'translate(' + horz + ',' + vert + ')';
     });
      tmp_leg.append('rect')
      .attr('width', legendRectSize)
      .attr('height', legendRectSize)
      .style('fill', function (d) {
       var col = colorScale(d);
       return d3.rgb(col).brighter([0.35]).toString();
     })
      .style('stroke', function (d) {
       var col = colorScale(d);
       return d3.rgb(col).toString();
     });

      tmp_leg.append('text')
      .attr('x', legendRectSize + legendSpacing + 5)
      .attr('y', legendRectSize - legendSpacing)
      .text(function(d) {
       return d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');
     });
  }
  return tmp_leg;
}

function removeLegend(){
	legend.selectAll('.legend')
  .remove();
    legend.selectAll('.map_feature')
  .remove();
}
var map_features;
function make_map(){
    var width = 300,
        height = 250,
        active = d3.select(null);

    console.log('enter map');
    var projection = d3.geo.robinson()
            .scale(3000)
            .translate([width, height * 1.5])

    var path = d3.geo.path()
        .projection(projection);

    var svg = d3.select("#legend");

    svg
        .attr("width", width)
        .attr("height", height);

    svg.append("rect")
        .attr("class", "map_background")
        .attr("width", width)
        .attr("height", height);

    var g = svg.append("g")
        .style("stroke-width", "1.5px");

    d3.json("/data/ebov.json", function(error, locations) {
        var locationData = topojson.feature(locations, locations.objects.ebov).features;

        map_features = g.selectAll(".map_feature")
            .data(locationData)
            .enter().append("path")
            .style("fill", function(d) {console.log(d.id);
                return (d.properties.ISO === "GIN" ? "lightseagreen" :
                        (d.properties.ISO === "SLE" ? "steelblue" :
                        (d.properties.ISO === "LBR" ? "lightcoral" : "lightseagreen")));
                })
            .attr("d", path)
            .attr("class", "map_feature")
            .on("mouseover",mouseOverMap)
            .on("mouseout",mouseOutMap);

        g.attr("transform", "translate(250,200) scale(" + 0.5 + ")");

      g.append("path")
          .datum(topojson.mesh(locations, locations.objects.ebov, function(a, b) { return a !== b; }))
          .attr("class", "map_mesh")
          .attr("d", path);

    });
}


function mouseOverMap(region){
    treeplot.selectAll(".tip")
            .filter(function (d){return d.region==region.id;})
                .attr("r", function(d){return tipRadius*2;});
}

function mouseOutMap(region){
    treeplot.selectAll(".tip")
            .filter(function (d){return d.region==region.id;})
                .attr("r", function(d){return tipRadius;});
}
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

function patch_region_name(d){
  var tmp = d.properties.NAME_2;
  if (tmp == null){
    tmp= d.properties.NAME_1;
  }
  if (tmp=='?'||tmp==null){
    tmp = d.properties.ISO;
    console.log('Falling back on ISO: '+ d.properties.NAME_2+ ' ' + d.properties.NAME_1 + ' ' + d.properties.ISO + ' ' +d.id);
  }
  return tmp.replace(' ','');
}

function patch_color(d) {
	return regionColorScale(patch_region_name(d));
}

function make_map(){
    var width = 400,
        height = 340,
        active = d3.select(null);

    console.log('enter map');
    var projection = d3.geo.robinson()
            .scale(3000)
            .translate([width, height * 1.5])

    var path = d3.geo.path()
        .projection(projection);

    var svg = d3.select("#legend");
    svg.call(mapTooltip);

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
            .style("fill", patch_color)
            .attr("d", path)
            .attr("class", "map_feature")
            .on("mouseover",mouseOverMap)
            .on("mouseout",mouseOutMap);

        g.attr("transform", "translate(260,180) scale(" + 0.7 + ")");

      g.append("path")
          .datum(topojson.mesh(locations, locations.objects.ebov, function(a, b) { return a !== b; }))
          .attr("class", "map_mesh")
          .attr("d", path);

    });
}

function match_region(map_region, tip){
  var tmp = patch_region_name(map_region);
  return tmp==tip.region.replace(' ','');
}

function patch_in_list(map_region, branch){
  var tmp = patch_region_name(map_region);
  for (reg in branch.target.allRegions){
    if (tmp==reg.replace(' ','')) return true;
  }
  return false;
}


function mouseOverMap(region){
	console.log(region);
    mapTooltip.show(region);
    treeplot.selectAll(".tip")
            .filter(function (d){ return match_region(region, d);})
                .attr("r", function(d){return tipRadius*2;})
                .style("fill", function (t) {
                  return d3.rgb(tipFillColor(t)).brighter();
                });
	legend.selectAll('.map_feature')
		.filter(function (m) { return patch_region_name(m) == patch_region_name(region);})
		.style("fill", function(m) {
			return d3.rgb(patch_color(region)).brighter();
		});
}

function mouseOutMap(region){
    mapTooltip.hide(region);
    treeplot.selectAll(".tip")
            .filter(function (d){ return match_region(region, d);})
                .attr("r", function(d){return tipRadius;})
                .style("fill", tipFillColor);
	legend.selectAll('.map_feature')
		.filter(function (m) { return patch_region_name(m) == patch_region_name(region);})
		.style("fill", function(m) {
			return d3.rgb(patch_color(region));
		});          
}
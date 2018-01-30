var virusTooltip = d3.tip()
	.direction('se')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {

		string = "";

		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";

		string += "<div class=\"smallnote\">";

		// check if vaccine strain
		if (vaccineStrains.indexOf(d.strain) != -1) {
			string += "Vaccine strain<br>";
			var vaccine_date = new Date(vaccineChoice[d.strain]);

			string += "First chosen " + vaccine_date.toLocaleString("en-us", { month: "short" }) + " " + vaccine_date.getFullYear() + "<br>";
			string += "<div class=\"smallspacer\"></div>";
		}

		// division / country / region
		// known to division level
		if (typeof d.attr.division != "undefined" && typeof d.attr.country != "undefined" && d.attr.division != d.attr.country) {
			string += d.attr.division.replace(/([A-Z])/g, ' $1').replace(/_/g, ' ').toTitleCase();
			string += ", ";
			string += d.attr.country.replace(/([A-Z])/g, ' $1').replace(/_/g, ' ').toTitleCase().replace(/^Usa/, 'USA').replace(/^Usvi/, 'USVI');
		}

		// division / country / region
		// known to country level
		else if (typeof d.attr.division != "undefined" && typeof d.attr.country != "undefined" && d.attr.division == d.attr.country) {
			string += d.attr.country.replace(/([A-Z])/g, ' $1').replace(/_/g, ' ').toTitleCase().replace(/^Usa/, 'USA').replace(/^Usvi/, 'USVI');
		}

		// country / region
		// known to country level
		else if (typeof d.attr.country != "undefined") {
			string += d.attr.country.replace(/([A-Z])/g, ' $1').replace(/_/g, ' ').toTitleCase().replace(/^Usa/, 'USA');
		}

		if (typeof d.attr.raw_date != "undefined") {
			string += ", " + d.attr.raw_date;
		}
		else if (typeof d.attr.date != "undefined") {
			string += ", " + d.attr.date;
		}
		if ((typeof d.attr.db != "undefined") && (typeof d.attr.accession != "undefined") && (d.attr.db == "GISAID")) {
			string += "<br>GISAID ID: EPI" + d.accession;
		}
		if ((typeof d.attr.db != "undefined") && (typeof d.attr.accession != "undefined") && (d.attr.db == "Genbank")) {
			string += "<br>Accession: " + d.accession;
		}
		if (typeof d.attr.lab != "undefined") {
			if (d.attr.lab != "") {
				string += "<br>Source: " + d.attr.lab.replace(/([A-Z])/g, ' $1').replace(/_/g, ' ').toTitleCase().substring(0,25);
				if (d.attr.lab.length>25) string += '...';
			}
		}
		if (typeof d.attr.authors != "undefined") {
			if ((d.attr.authors != "") && (d.attr.authors != "?")) {
				string += "<br>Authors: " + d.attr.authors.substring(0,25);
				if (d.attr.authors.length>25) string += '...';
			}
		}
		string += "</div>";
		// following may or may not be present
		if ((typeof focusNode != "undefined")){
			string += "<div class=\"smallspacer\"></div>";
			string += "HI against serum from "+focusNode.strain;
			string += "<div class=\"smallspacer\"></div>";
			string += "<div class=\"smallnote\">"
			string += '<table class="table table-condensed"><thead><tr><td>Serum</td><td>&#916log<sub>2</sub></td><td>heterol.</td><td>homol.</td></tr></thead><tbody>';
			var tmp_titers = HI_titers[focusNode.clade][d.clade];
			var tmp_auto_titers = HI_titers[focusNode.clade][focusNode.clade];
			var tmp_avi = titer_subs_model["avidity"][d.clade];
			var tmp_pot = titer_subs_model["avidity"][d.clade];
			if (typeof tmp_titers != "undefined"){
				for (var tmp_serum in tmp_titers){
					var autoHI = tmp_auto_titers[tmp_serum][1];
					var rawHI = tmp_titers[tmp_serum][1];
					var logHI = tmp_titers[tmp_serum][0];
					if (correctVirus){logHI-=tmp_avi;}
					if (correctPotency){logHI-=titer_subs_model["potency"][focusNode.clade][tmp_serum];}
					var serum_name;
					if (tmp_serum.length<20){
						serum_name = tmp_serum;
					}else{
						serum_name = tmp_serum.substring(0,17)+'...';
					}
					string += '<tr><td>' + serum_name + '</td><td>' +  logHI.toFixed(2) + '</td><td>' + rawHI.toFixed(0)+ '</td><td>' + autoHI.toFixed(0) +"</td></tr>";
				}
			}
			string += '<tr><td>' + 'Tree model' + '</td><td>' +  d.HI_dist_tree.toFixed(2) + '</td><td> --- </td><td>---</td></tr>';
			string += '<tr><td>' + 'Subs. model ' + '</td><td>' +  d.HI_dist_mut.toFixed(2) + '</td><td> --- </td><td>---</td></tr>';
			string += "</tbody></table></div>";
		}

		string += "<div class=\"smallspacer\"></div>";
		// following may or may not be present
		string += "<div class=\"smallnote\">";
		if (typeof d.attr.cTiter != "undefined") {
			string += "Antigenic adv: " + d.attr.cTiter.toFixed(1) + "<br>";
		}
		if (typeof d.attr.ep != "undefined") {
			string += "Epitope distance: " + d.attr.ep + "<br>";
		}
		if (typeof d.attr.rb != "undefined") {
			string += "Receptor binding distance: " + d.attr.rb + "<br>";
		}
		if (typeof d.attr.glyc != "undefined") {
			string += "Glycosylation sites: " + d.attr.glyc + "<br>";
		}
		if (typeof d.attr.age != "undefined") {
			if (d.attr.age != "unknown"){
				string += "Host age: " + d.attr.age + "y<br>";
			}
		}
		if (typeof d.LBI != "undefined") {
			string += "Local branching index: " + d.LBI.toFixed(3) + "<br>";
		}
		if (typeof d.dfreq != "undefined") {
			string += "Freq. change: " + d.dfreq.toFixed(3) + "<br>";
		}
		if (typeof d.attr.pred_distance != "undefined") {
			string += "Predicted distance: " + d.pred_distance.toFixed(3) + "<br>";
		}
		string += "</div>";
		return string;
	});
treeplot.call(virusTooltip);


var linkTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
		}
		if (typeof d.attr.dTiter != "undefined") {
			string += "<br>Titer drop: " + d.attr.dTiter.toFixed(2)
		}
		string += "<div class=\"smallspacer\"></div>";
		string += "<div class=\"smallnote\">";
		var ncount = d.muts.length;
		if (ncount) {string += "<b>Mutations:</b><ul>";}
		if ((typeof d.aa_muts !="undefined")){
			for (tmp_gene in d.aa_muts) {ncount+=d.aa_muts[tmp_gene].length;}
			for (tmp_gene in d.aa_muts){
				if (d.aa_muts[tmp_gene].length){
					string+="<li>"+tmp_gene+":</b> "+d.aa_muts[tmp_gene].join(', ') + "</li>";
				}
			}
		}
		if ((typeof d.muts !="undefined") && (d.muts.length)){
			var nmuts = d.muts.length;
			var tmp_muts = d.muts.slice(0,Math.min(10, nmuts))
			string += "<li>nuc: "+tmp_muts.join(', ');
			if (nmuts>10) {string+=' + '+ (nmuts-10) + ' more';}
			string += "</li>";
		}
		string += "</ul>";
		if (typeof d.fitness != "undefined") {
			string += "Fitness: " + d.fitness.toFixed(3) + "<br>";
		}
		if (typeof d.attr.age_score != "undefined") {
			string += "Age_stat: " + d.attr.age_score.toFixed(3) + "<br>";
		}
		if (typeof d.attr.avg_age != "undefined") {
			string += "Avg age: " + d.attr.avg_age.toFixed(3) + "y<br>";
		}

		string += "click to zoom into clade"
		string += "</div>";
		return string;
	});
treeplot.call(linkTooltip);


var matchTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = d.strain+ "<i> is closest match of:</i><ul>";
		string += "<div class=\"smallspacer\"></div>";
		for (var mi=0; mi<d.matches.length;mi++){
			string+="<li>" +d.matches[mi].substring(0,Math.min(30,d.matches[mi].length))+'</li>';
		}
		string += "</ul>";
		return string;
	});
treeplot.call(matchTooltip);

var mapTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0,10])
	.html(function(d) {
		  var tmp = '';
		  if (d.properties.NAME_2!=null) tmp+=d.properties.NAME_2;
		  if ((d.properties.NAME_2==null)&&(d.properties.NAME_1!=null)){
		  	tmp+=d.properties.NAME_1;
		  }
		  tmp+=', '+d.properties.ISO;
			tmp = tmp.replace(/_/g, ' ').toTitleCase();
		  return tmp;
	});

// configure mathjax to work with markdown through script tags: `

MathJax.Hub.Config({
	tex2jax: { skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'], inlineMath: [["$","$"],["\\(","\\)"]] }
});
	
MathJax.Hub.Queue(function() {
	var all = MathJax.Hub.getAllJax(), i;
	for(i=0; i < all.length; i += 1) {
		all[i].SourceElement().parentNode.className += ' has-jax';
	}
});

// Switch MathJax font to sans serif

MathJax.Hub.Config({
  "HTML-CSS": {availableFonts: ["TeX"]},
	  MMLorHTML: {prefer: "HTML"}
});

MathJax.Hub.Register.StartupHook("HTML-CSS Jax Ready",function () {
  var VARIANT = MathJax.OutputJax["HTML-CSS"].FONTDATA.VARIANT;
  VARIANT["normal"].fonts.unshift("MathJax_SansSerif");
  VARIANT["bold"].fonts.unshift("MathJax_SansSerif-bold");
  VARIANT["italic"].fonts.unshift("MathJax_SansSerif-italic");
  VARIANT["-tex-mathit"].fonts.unshift("MathJax_SansSerif-italic");
});

MathJax.Hub.Register.StartupHook("SVG Jax Ready",function () {
  var VARIANT = MathJax.OutputJax.SVG.FONTDATA.VARIANT;
  VARIANT["normal"].fonts.unshift("MathJax_SansSerif");
  VARIANT["bold"].fonts.unshift("MathJax_SansSerif-bold");
  VARIANT["italic"].fonts.unshift("MathJax_SansSerif-italic");
  VARIANT["-tex-mathit"].fonts.unshift("MathJax_SansSerif-italic");
});
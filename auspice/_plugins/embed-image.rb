# Author: Trevor Bedford
# License: MIT

# Examples:
#  {% embed_image url %}

module Jekyll
	class EmbedImage < Liquid::Tag
		def initialize(tag_name, markup, tokens)
			super
			@markup = "#{markup}".strip
		end
		def render(context)
		
			parsed = Liquid::Template.parse(@markup).render context
			url = parsed.split(/ /).first
			url.gsub!(/ /, '%20')	
			html = ""
			html += "<div class=\"row\">"
			html += "<div class=\"col-lg-1\"></div>"
			html += "<div class=\"col-lg-10\">"
			html += "<img src=\"#{url}\" class=\"img-responsive\"/>"
			html += "</div>"
			html += "<div class=\"col-lg-1\"></div>"
			html += "</div>"
			html 
			
		end
	end
end

Liquid::Template.register_tag('embed_image', Jekyll::EmbedImage)

# Author: Trevor Bedford
# License: MIT

# Examples:
#  {% embed_image url width %}
# width should be an even number between 2 and 12 for the number of columns the image takes up

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
			width = parsed.split(/ /).drop(1).join(' ')
			width.gsub!(/ /, '%20')
			width = width.to_i
			border = (12-width)/2
			
			html = ""
			html += "<div class=\"spacer\"></div>"
			html += "<div class=\"row\">"
			if border > 0 then
				html += "<div class=\"col-lg-#{border}\"></div>"
			end
			html += "<div class=\"col-lg-#{width}\">"
			html += "<a href=\"#{url}\" target=\"_blank\">"
			html += "<img src=\"#{url}\" class=\"img-responsive\"/>"
			html += "</a>"
			html += "</div>"
			if border > 0 then
				html += "<div class=\"col-lg-#{border}\"></div>"
			end
			html += "</div>"
			html += "<div class=\"spacer\"></div>"			
			html 
			
		end
	end
end

Liquid::Template.register_tag('embed_image', Jekyll::EmbedImage)

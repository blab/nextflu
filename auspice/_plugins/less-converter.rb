# less converter via https://gist.github.com/KBalderson/5689220

module Jekyll
  class LessConverter < Converter
    safe true
    priority :high
    
    def setup
      return if @setup
      require 'less'
      @setup = true
    rescue LoadError
      STDERR.puts 'You are missing the library required for less. Please run:'
      STDERR.puts ' $ [sudo] gem install less'
      raise FatalException.new("Missing dependency: less and/or therubyracer")
    end
    
    def matches(ext)
      ext =~ /less|lcss/i
    end
    
    def output_ext(ext)
      ".css"
    end
    
	def convert(content)
      setup
      begin
        parser = Less::Parser.new
        parser = parser.parse(content).to_css
      rescue => e
        puts "Less Exception: #{e.message}"
      end
    end

  end
end
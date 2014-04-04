# network visualization of the CTC 2013 abstracts

# the function that does all of the work
draw = (data) -> 
  # Width and height of SVG
  w = 750
  h = 750
  pad = {left:60, top:40, right:40, bottom: 40}

  darkBlue = "darkslateblue"
  green = "#060"
  gray = "#aaa"
  red = "crimson"

  # source/target version of edges
  data.edges2 = []
  for e in data.edges
    data.edges2.push({source:e[0],target:e[1]})

  # min and max nauthors
  minAuthors = 50
  maxAuthors = 0
  for a in data.abstracts
    minAuthors = a.nauthors if a.nauthors < minAuthors
    maxAuthors = a.nauthors if a.nauthors > maxAuthors

  # tool tip
  tip = d3.svg.tip()
          .orient("right")
          .padding(3)
          .text((d) -> "#{d.type} #{d.number}")
          .attr("class", "d3-tip")
          .attr("id", "tip")

  # Initialize a default force layout, using the nodes and edges in data
  force = d3.layout.force()
                   .nodes(data.abstracts)
                   .links(data.edges2)
                   .size([w, h])
                   .linkDistance([50])
                   .charge([-100])
                   .start()

  colors = d3.scale.category10()

  svg = d3.select("div#abstracts")
          .append("svg")
          .attr("width", w+pad.left+pad.right)
          .attr("height", h+pad.top+pad.bottom)

  svg.append("rect")
     .attr("x", pad.left)
     .attr("y", pad.top)
     .attr("width", w)
     .attr("height", h)
     .attr("fill", "none")
     .attr("stroke", "black")
     .attr("stroke-width", 2)
      
  # color points by number of abstracts
  colorScale = d3.scale.log()
                 .domain([1, maxAuthors])
                 .range(["white", darkBlue])
                 .clamp(true)

  # Create edges as lines
  edges = svg.selectAll("line")
             .data(data.edges2)
             .enter()
             .append("line")
             .style("stroke", gray)
             .style("stroke-width", 1)
      
  # Create nodes as circles
  nodes = svg.selectAll("circle")
             .data(data.abstracts)
             .enter()
             .append("a")
             .attr("xlink:href", (d) -> "#{d.file}")
             .append("circle")
             .attr("r", 10)
             .style("fill", (d) -> colorScale(d.nauthors))
             .style("stroke", "black")
             .style("stroke-width", 1)
             .call(force.drag)
             .on "mouseover", (d) ->
                 d3.select(this).attr("r",15)
                 d3.select("body").selectAll("p#title").text(d.title)
                 tip.call(this, d)
             .on "mouseout", (d) ->
                 d3.select(this).attr("r",10)
                 d3.select("body").selectAll("p#title").text("")
                 d3.selectAll("#tip").remove()
             .on "click", (d) ->
                 d3.select(this).attr("r",10)
                 d3.select("body").select("p#title").text("")
                 d3.selectAll("#tip").remove()

  # Every time the simulation "ticks", this will be called
  force.on("tick", ->
              edges.attr("x1", (d) -> d.source.x+pad.left)
                   .attr("y1", (d) -> d.source.y+pad.top)
                   .attr("x2", (d) -> d.target.x+pad.left)
                   .attr("y2", (d) -> d.target.y+pad.top)
              nodes.attr("cx", (d) -> d.x+pad.left)
                   .attr("cy", (d) -> d.y+pad.top))

# load json file and call draw function
d3.json("abstracts.json", draw)

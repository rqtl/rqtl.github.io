# network visualization of the CTC 2013 abstracts
# here nodes = authors

# the function that does all of the work
draw = (data) -> 
  # Width and height of SVG
  w = 750
  h = 750
  pad = {left:60, top:40, right:40, bottom: 40}

  darkBlue = "darkslateblue"
  green = "#060"
  red = "crimson"
  gray = "#aaa"

  # source/target version of edges
  data.authoredges2 = []
  for e in data.authoredges
    data.authoredges2.push({source:e[0],target:e[1]})

  # min and max nabstracts
  minAbstracts = 50
  maxAbstracts = 0
  for a in data.authors
    minAbstracts = a.nabstracts if a.nabstracts < minAbstracts
    maxAbstracts = a.nabstracts if a.nabstracts > maxAbstracts

  # Initialize a default force layout, using the nodes and edges in data
  force = d3.layout.force()
                   .nodes(data.authors)
                   .links(data.authoredges2)
                   .size([w, h])
                   .linkDistance([20])
                   .charge([-50])
                   .start()

  svg = d3.select("div#authors")
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
      
  # tool tip
  tip = d3.svg.tip()
          .orient("right")
          .padding(3)
          .text((d) -> d.name)
          .attr("class", "d3-tip")
          .attr("id", "tip")

  # color points by number of abstracts
  colorScale = d3.scale.log()
                 .domain([1, maxAbstracts])
                 .range(["white", darkBlue])
                 .clamp(true)

  # Create edges as lines
  edges = svg.selectAll("line")
             .data(data.authoredges2)
             .enter()
             .append("line")
             .style("stroke", gray)
             .style("stroke-width", 1)
      
  # Create nodes as circles
  nodes = svg.selectAll("circle")
             .data(data.authors)
             .enter()
             .append("circle")
             .attr("r", 6)
             .style("fill", (d) -> colorScale(d.nabstracts))
             .style("stroke", "black")
             .style("stroke-width", 1)
             .call(force.drag)
             .on "mouseover", (d) ->
                  d3.select(this).attr("r", 8)
                  tip.call(this, d)
             .on "mouseout", (d) ->
                  d3.select(this).attr("r", 6)
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

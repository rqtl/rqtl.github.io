# programvis.coffee
# Simple interactive with program

# function that does the work
draw = (data) ->
  w = 798
  h = 850

  pad = {left:50, top:10, right:10, bottom:20, rect:50}
  tickWidth = 4
  totalw = w + pad.left + pad.right
  totalh = h + pad.top + pad.bottom

  vlightGray = "#EEE"  
  lightGray = "#CCC"
  darkGray = "#AAA"
  pink = "#E9CFEC"
  purple = "#8C4374"
  borderColor = "black"
  hlineColor = "black"
  textColor = "black"
  borderWidth = 1
  thickBorderWidth = 4
  borderHilit = d3.rgb(0, 0, 128)
  hlineWidth = 1
  headHeight = 50
  scienceColor = d3.rgb(201,255,255)
  keynoteColor = d3.rgb(121,255,255)
  partyColor = d3.rgb(255,201,237)

  nDays = data.days.length
  rectw = (w - pad.rect*(nDays-1))/nDays
  rectLeft = (i*(rectw+pad.rect) for i of data.days)

  svg = d3.select("div#figure").append("svg")
          .attr("width", totalw)
          .attr("height", totalh)

  fig = svg.append("g").attr("id", "program")
           .attr("transform", "translate(#{pad.left},#{pad.top})")

  # find min and max times (as numbers)
  minTime = 24
  maxTime = 0
  for activity in data.activities
    minTime = activity.startnum if activity.startnum < minTime
    maxTime = activity.endnum if activity.endnum > maxTime

  yScale = d3.scale.linear()
                   .domain([minTime, maxTime])
                   .range([headHeight, h])

  # column rectangles
  fig.selectAll("empty")
     .data(data.days)
     .enter()
     .append("a")
     .attr("xlink:href", (d,i) -> "program_#{data.sdays[i]}.html")
       .append("rect")
       .attr("x", (d,i) -> rectLeft[i])
       .attr("y", 0)
       .attr("height", h)
       .attr("width", rectw)
       .attr("stroke", borderColor)
       .attr("fill", vlightGray)
       .attr("stroke-width", borderWidth)
       .on "mouseover", (d,i) ->
           d3.selectAll("rect#border#{i}")
             .attr("stroke-width", thickBorderWidth)
             .attr("stroke", borderHilit)
       .on "mouseout", (d,i) ->
           d3.selectAll("rect#border#{i}")
             .attr("stroke-width", borderWidth)
             .attr("stroke", borderColor)
       .on "click", (d,i) ->
           d3.selectAll("rect#border#{i}")
             .attr("stroke-width", borderWidth)
             .attr("stroke", borderColor)

  # head rectangles
  fig.selectAll("empty")
     .data(data.days)
     .enter()
     .append("rect")
     .attr("x", (d,i) -> rectLeft[i])
     .attr("y", 0)
     .attr("height", headHeight)
     .attr("width", rectw)
     .attr("stroke", borderColor)
     .attr("fill", darkGray)
     .attr("stroke-width", borderWidth)
     .attr("id", (d,i) -> "headrect#{i}")
     .style("pointer-events", "none")

  # head text
  fig.selectAll("empty")
     .data(data.days)
     .enter()
     .append("text")
     .text((d) -> d)
     .attr("x", (d,i) -> rectLeft[i] + rectw/2)
     .attr("y", headHeight/2)
     .attr("fill", textColor)
     .attr("dominant-baseline", "middle")
     .attr("text-anchor", "middle")
     .attr("id", (d,i) -> "headtext#{i}")
     .style("font-size", "0.8em")
     .style("pointer-events", "none")

  # activity rectangles
  fig.selectAll("empty")
     .data(data.activities)
     .enter()
     .append("rect")
     .attr("x", (d) -> rectLeft[d.day])
     .attr("width", rectw)
     .attr("y", (d) -> yScale(d.startnum))
     .attr("height", (d) -> yScale(d.endnum) - yScale(d.startnum))
     .attr("stroke", borderColor)
     .attr("stroke-width", borderWidth)
     .attr("fill", (d) ->
              return scienceColor if d.type == "science"
              return keynoteColor if d.type == "keynote"
              return partyColor if d.type == "party"
              lightGray)
     .style("pointer-events", "none")

  # activity text
  fig.selectAll("empty")
     .data(data.activities)
     .enter()
     .append("text")
     .text((d) -> d.activity)
     .attr("x", (d) -> rectLeft[d.day] + rectw/2)
     .attr("y", (d) -> (yScale(d.endnum) + yScale(d.startnum))/2)
     .attr("fill", textColor)
     .attr("dominant-baseline", "middle")
     .attr("text-anchor", "middle")
     .style("font-size", "0.75em")
     .style("pointer-events", "none")

  # change 12:00 to noon
  for a in data.activities
    a.start = "noon" if a.start == "12:00"
    a.end = "noon" if a.end == "12:00"

  # start labels
  fig.selectAll("empty")
     .data(data.activities)
     .enter()
     .append("text")
     .text((d) -> d.start)
     .attr("y", (d) -> yScale(d.startnum))
     .attr("x", (d) -> rectLeft[d.day] - tickWidth*2)
     .attr("fill", textColor)
     .attr("dominant-baseline", "middle")
     .attr("text-anchor", "end")
     .style("font-size", "0.7em")
     .style("pointer-events", "none")

  # start ticks
  fig.selectAll("empty")
     .data(data.activities)
     .enter()
     .append("line")
     .attr("y1", (d) -> yScale(d.startnum))
     .attr("y2", (d) -> yScale(d.startnum))
     .attr("x1", (d) -> rectLeft[d.day] - tickWidth)
     .attr("x2", (d) -> rectLeft[d.day])
     .attr("stroke", hlineColor)
     .attr("stroke-width", hlineWidth)
     .style("pointer-events", "none")

  needs_endlabel = []
  for a in data.activities
    needs_endlabel.push(a) if a.activity == "Reception" or a.activity == "Social night out" or
                              a.activity == "Session 5" or a.activity == "Posters (odd #\u2019s)"

  # end labels
  fig.selectAll("empty")
     .data(needs_endlabel)
     .enter()
     .append("text")
     .text((d) -> d.end)
     .attr("y", (d) -> yScale(d.endnum))
     .attr("x", (d) -> rectLeft[d.day] - tickWidth*2)
     .attr("fill", textColor)
     .attr("dominant-baseline", "middle")
     .attr("text-anchor", "end")
     .style("font-size", "0.7em")
     .style("pointer-events", "none")

  # end ticks
  fig.selectAll("empty")
     .data(needs_endlabel)
     .enter()
     .append("line")
     .attr("y1", (d) -> yScale(d.endnum))
     .attr("y2", (d) -> yScale(d.endnum))
     .attr("x1", (d) -> rectLeft[d.day] - tickWidth)
     .attr("x2", (d) -> rectLeft[d.day])
     .attr("stroke", hlineColor)
     .attr("stroke-width", hlineWidth)
     .style("pointer-events", "none")

  # borders again
  fig.selectAll("empty")
     .data(data.days)
     .enter()
     .append("rect")
     .attr("x", (d,i) -> rectLeft[i])
     .attr("y", 0)
     .attr("height", h)
     .attr("width", rectw)
     .attr("stroke", borderColor)
     .attr("fill", "none")
     .attr("stroke-width", borderWidth)
     .attr("id", (d,i) -> "border#{i}")
     .style("pointer-events", "none")

# load json file and call draw function
d3.json("programvis/schedule.json", draw)

---
layout: archive
title: "Publications"
permalink: /publications/
author_profile: true
---

{% include base_path %}

{% if "a" == "a" %}
  You can also find my articles on <u><a href="https://scholar.google.com/citations?user=IJMEC2EAAAAJ&hl=en">my Google Scholar profile</a>.</u>
{% endif %}

<p><u>S. Sun</u>, <u>I. Eisenman</u>, and A. Stewart (2016). 
<b>The influence of Southern Ocean surface buoyancy forcing on glacial-interglacial changes in the global deep ocean stratification.</b> 
<i>Geophys Res Lett</i> 43, 8124-8132.


{% for post in site.publications reversed %}
  {% include archive-single.html %}
{% endfor %}

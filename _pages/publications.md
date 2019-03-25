---
layout: archive
title: "Publications"
permalink: /publications/
author_profile: true
---

{% include base_path %}

{% if "a" == "a" %}
   my google scholar: {{ author.name }}
{% endif %}

{% if "a" == "a" %}
  You can also find my articles on <u><a href="https://scholar.google.com/citations?user=IJMEC2EAAAAJ&hl=en">my Google Scholar profile</a>.</u>
{% endif %}

{% for post in site.publications reversed %}
  {% include archive-single.html %}
{% endfor %}

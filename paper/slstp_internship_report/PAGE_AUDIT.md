# Page and rendering audit

Final artifact: `manuscript.pdf`

Handoff rerender checked: August 13, 2026 UTC.

- Paper size: A4
- Total pages: 11
- Scientific main text: pages 1-9
- References begin: page 9, after the main text and report back matter
- Appendix begins: page 10, after the references
- Main-text limit: passed; the scientific narrative is under 15 pages before references and appendix

All 11 pages were rendered to PNG and inspected at full-page scale after the final edit. Dense figure and table pages were also checked individually at full resolution. Page 1 was checked after shortening the abstract. Pages 2 through 6 were checked for the GLARE, expiMap, and generator figures. Pages 7 through 9 were checked after condensing the tissue hypotheses and discussion. Pages 9 and 10 were checked at the transition from the main text to the references and appendix. Page 11 was checked after the appendix was made self-contained. The audit found:

- no clipped or overlapping text;
- no split tables;
- no missing or blank figures;
- readable axes, labels, legends, and captions;
- consistent heading, body, caption, and table styles;
- 9 main figures and 2 appendix figures rendered in their intended order;
- no malformed replacement characters in extracted PDF text;
- all 45 numbered references cited in the main text.

Rebuild command:

```bash
python \
  -m nasa_mouse_internship_report.build_report
```

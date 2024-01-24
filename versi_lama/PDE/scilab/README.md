Convert `*.png` to `.gif`:

```bash
convert -delay 10 -loop 0 TEMP_*.png anim01.gif
```

Convert `pdf` to `png`, with resize

```bash
convert -density 300 a.pdf -resize 25% a.png
```

Try using `-units PixelsPerInch`


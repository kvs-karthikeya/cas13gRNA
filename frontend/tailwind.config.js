/** Tailwind build for the Cas13 gRNA Designer front end.
 *  Rebuild after editing index.html:   npx tailwindcss -i src/input.css -o assets/tailwind.css --minify */
module.exports = {
  darkMode: 'class',
  content: ['./index.html'],
  theme: {
    extend: {
      fontFamily: {
        sans: ['"Open Sans"', 'ui-sans-serif', 'system-ui', 'sans-serif'],
        mono: ['"JetBrains Mono"', 'ui-monospace', 'SFMono-Regular', 'Menlo', 'monospace'],
      },
    },
  },
  plugins: [],
};

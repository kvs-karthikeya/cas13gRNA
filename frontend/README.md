# Cas13 gRNA Designer - front end

Static site: `index.html` + `assets/`. Nothing to build to deploy it.

1. Open `index.html` and set the backend address in the first script tag:
   `window.CAS13_API = 'https://YOUR-APP.up.railway.app';`   (https, no trailing slash)
2. Deploy this folder (Netlify: publish directory = `frontend`).

Black and white themes are switched with the sun/moon button; the choice is remembered in the browser.

Only if you edit the HTML classes: `npm install` then `npm run build:css` to regenerate `assets/tailwind.css`.
Colors, glow and component styles live in `assets/theme.css`.

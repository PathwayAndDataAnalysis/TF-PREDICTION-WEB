// tf-pred-webserver/tailwind.config.js
module.exports = {
    content: [
        "./templates/**/*.html",
        "./app/**/*.py",
        "./static/css/src/**/*.css",
        "./static/js/**/*.js",
    ],
    theme: {
        extend: {
            colors: {
                blue: {
                  50:  '#f2f4f4',
                  100: '#e1e5e5',
                  200: '#c3cbcb',
                  300: '#a4b0b1',
                  400: '#859596',
                  500: '#4a5759', // base
                  600: '#3f4b4c',
                  700: '#343f40',
                  800: '#2a3334',
                  900: '#1f2728',
                },
            },
        },
    },
    plugins: [],
}

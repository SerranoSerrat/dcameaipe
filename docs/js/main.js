// ---------- Mobile sidebar toggle ----------
const navToggle = document.querySelector('.nav-toggle');
const scrim = document.querySelector('.sidebar-scrim');

function closeNav() { document.body.classList.remove('nav-open'); }
function toggleNav() { document.body.classList.toggle('nav-open'); }

if (navToggle) navToggle.addEventListener('click', toggleNav);
if (scrim) scrim.addEventListener('click', closeNav);

// ---------- Highlight the current page in the sidebar ----------
(function highlightActive() {
  const here = location.pathname.split('/').pop() || 'index.html';
  document.querySelectorAll('.sidebar a[href]').forEach((link) => {
    const href = link.getAttribute('href').split('/').pop();
    if (href === here) link.classList.add('active');
  });
})();

// ---------- Sidebar search filter ----------
const searchInput = document.querySelector('.sidebar-search input');
if (searchInput) {
  searchInput.addEventListener('input', (e) => {
    const q = e.target.value.trim().toLowerCase();
    document.querySelectorAll('.sidebar li').forEach((li) => {
      const text = li.textContent.toLowerCase();
      li.style.display = !q || text.includes(q) ? '' : 'none';
    });
  });
}

// ---------- Copy buttons on code blocks ----------
document.querySelectorAll('.code-block').forEach((block) => {
  const btn = document.createElement('button');
  btn.className = 'copy-btn';
  btn.type = 'button';
  btn.textContent = 'Copy';
  btn.addEventListener('click', async () => {
    const code = block.querySelector('code');
    try {
      await navigator.clipboard.writeText(code.innerText);
      btn.textContent = 'Copied';
      setTimeout(() => { btn.textContent = 'Copy'; }, 1500);
    } catch (err) {
      btn.textContent = 'Select & copy';
    }
  });
  block.appendChild(btn);
});
